import logging
import os
import uuid
from collections import OrderedDict, Counter

from kb_reaction_gene_finder.core.re_api import RE_API
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.KBaseReportClient import KBaseReport


class AppImpl:
    def __init__(self, config, ctx):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.re_api = RE_API(config['re-api'], ctx['token'])
        self.gfu = GenomeFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)

    @staticmethod
    def _validate_params(params, required, optional=set()):
        """Validates that required parameters are present. Warns if unexpected parameters appear"""
        required = set(required)
        optional = set(optional)
        pkeys = set(params)
        if required - pkeys:
            raise ValueError("Required keys {} not in supplied parameters"
                             .format(", ".join(required - pkeys)))
        defined_param = required | optional
        for param in params:
            if param not in defined_param:
                logging.warning("Unexpected parameter {} supplied".format(param))

    @staticmethod
    def _make_fasta(sequences, file_path):
        """Make a FASTA file from RE query output"""
        if not sequences or not isinstance(sequences, list) or "key" not in sequences[0] \
                or "seq" not in sequences[0]:
            raise ValueError(f"Unexpected input type:{sequences}")
        with open(file_path, 'w') as outfile:
            for seq in sequences:
                if seq["seq"]:
                    outfile.write(f'>{seq["key"]}\n{seq["seq"]}\n')

    def _find_best_homologs(self, query_seq_file, target_seq_file, rxn_id, noise_level=50,
                            number_vals_to_report=5, threads=1):
        """Blast the query_seq_file against the target_seq_file and return the best hits"""
        logging.info("running blastp for {0} vs {1}".format(query_seq_file, target_seq_file))
        out_cols = ['sseqid', 'qseqid', 'bitscore', 'pident', 'length', 'mismatch', 'evalue',]
        col_names = ['Genome Gene', 'Closest Database Gene', 'Bit Score', 'Percent Identity',
                     'Match Length', 'Mismatches', 'E Value']
        # wasn't able to get a pipe going
        # proc = subprocess.run( 'blastp -outfmt 6 -subject ' + target_seq_file +' -query ' + query_seq_file,
        #                       stdout=subprocess.PIPE, shell=True, universal_newline=True )
        #
        # for line in proc.stdout:
        #    print( line )

        # so using file output instead
        tmp_blast_output_file = os.path.join(self.scratch, "blastp.results" + str(uuid.uuid4()))

        blastp_cmd = f'blastp -outfmt "6 {" ".join(out_cols)}" -subject {target_seq_file} '\
                     f'-num_threads {threads} -query {query_seq_file} > {tmp_blast_output_file}'
        os.system(blastp_cmd)

        gene_hits = dict()
        top_bitscore = Counter()
        gene_hit_count = Counter()

        with open(tmp_blast_output_file) as bl:
            for line in bl:
                cols = OrderedDict(zip(col_names, line.strip().split()))
                bl_score = float(cols['Bit Score'])
                if bl_score < noise_level:
                    continue

                gene_hit_count[cols['Genome Gene']] += 1
                if cols['Genome Gene'] not in gene_hits or \
                        top_bitscore[cols['Genome Gene']] < bl_score:
                    gene_hits[cols['Genome Gene']] = cols
                    top_bitscore[cols['Genome Gene']] = float(cols['Bit Score'])

        top_genes = (gene[0] for gene in top_bitscore.most_common(number_vals_to_report))
        top_records = [{"Reaction ID": rxn_id, **gene_hits[gene],
                        "Total Gene Hits": str(gene_hit_count[gene])}
                       for gene in top_genes]

        return top_records

    def find_genes_from_similar_reactions(self, params):
        self._validate_params(params, {'workspace_name', 'query_genome_ref', },
                              {'number_of_hits_to_report', 'smarts_set', 'blast_score_floor',
                               'structural_similarity_floor', 'difference_similarity_floor',
                               'reaction_set', 'bulk_reaction_ids'})

        feature_seq_path = self.gfu.genome_proteins_to_fasta(
            {'genome_ref': params['query_genome_ref'],
             'include_functions': True,
             'include_aliases': False})['file_path']

        if not params.get('reaction_set'):
            params['reaction_set'] = []

        if params.get('bulk_reaction_ids'):
            params['reaction_set'] += params['bulk_reaction_ids'].split('\n')

        output = {'gene_hits': []}
        for rxn in params['reaction_set']:
            output['gene_hits'].extend(self.find_genes_for_rxn(rxn, feature_seq_path, params))
        output.update(self._build_report(params['reaction_set'],
                                         output['gene_hits'],
                                         params['workspace_name'],
                                         ))
        return output

    def find_genes_for_rxn(self, reaction, genome_feature_path, params):
        """Finds genes for a particular reaction using RE and BLAST"""
        search_seqs = self.re_api.get_related_sequences(
            reaction,
            params.get('structural_similarity_floor', 1),
            params.get('difference_similarity_floor', 1))
        if not search_seqs:
            return []

        search_fasta = f'rxn_search_{uuid.uuid4()}.fasta'
        self._make_fasta(search_seqs, search_fasta)

        return self._find_best_homologs(search_fasta,
                                        genome_feature_path,
                                        reaction,
                                        params.get('blast_score_floor', 50),
                                        params.get('number_of_hits_to_report', 5))

    def _build_report(self, reactions, gene_hits, workspace_name):
        """
                _generate_report: generate summary report for upload
                """
        output_html_files = self._generate_report_html(reactions, gene_hits)

        report_params = {
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            #'objects_created': [{'ref': pdb_obj_ref,
            #                     'description': 'Imported PDB'}],
            'workspace_name': workspace_name,
            'report_object_name': 'find_genes_for_rxn_' + str(uuid.uuid4())}

        output = self.kbr.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_report_html(self, reactions, gene_hits):
        """
            _generate_report: generates the HTML for the upload report
        """
        html_report = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'find_genes_for_rxn.html')

        # Build HTML tables for results
        table_lines = []
        table_lines.append(f'<h3 style="text-align: center">Best matching Genes</h3>')
        if not gene_hits:
            table_lines.append("<h5>No hits found!</h5>")
        else:
            table_lines.append('<table class="table table-bordered table-striped">')
            header = "</td><td>".join(gene_hits[0].keys())
            table_lines.append(f'\t<thead><tr><td>{header}</td></tr></thead>')
            table_lines.append('\t<tbody>')
            for row in gene_hits:
                line = "</td><td>".join(row.values())
                table_lines.append(f'\t\t<tr><td>{line}</td></tr>')
            table_lines.append('\t</tbody>')
            table_lines.append('</table>\n')

        # Fill in template HTML
        with open(os.path.join(os.path.dirname(__file__), 'find_genes_for_rxn_template.html')
                  ) as report_template_file:
            report_template = report_template_file.read() \
                .replace('*TABLES*', "\n".join(table_lines))

        with open(result_file_path, 'w') as result_file:
            result_file.write(report_template)

        html_report.append({'path': output_directory,
                            'name': os.path.basename(result_file_path),
                            'description': 'HTML report for reaction to gene candidates app'})

        return html_report
