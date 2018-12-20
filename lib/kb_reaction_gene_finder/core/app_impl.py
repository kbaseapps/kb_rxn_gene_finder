import logging
import os
import uuid
from collections import namedtuple

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

    @staticmethod
    def _find_best_homologs(query_seq_file, target_seq_file, noise_level=50,
                            number_vals_to_report=5, threads=1):
        """Blast the query_seq_file against the target_seq_file and return the best hits"""
        logging.info("running blastp for {0} vs {1}".format(query_seq_file, target_seq_file))
        blast_record = namedtuple('blast_record', ['qseqid', 'sseqid', 'pident', 'length',
                                                   'mismatch', 'gapopen', 'qstart', 'qend',
                                                   'sstart', 'send', 'evalue', 'bitscore'])
        # wasn't able to get a pipe going
        # proc = subprocess.run( 'blastp -outfmt 6 -subject ' + target_seq_file +' -query ' + query_seq_file,
        #                       stdout=subprocess.PIPE, shell=True, universal_newline=True )
        #
        # for line in proc.stdout:
        #    print( line )

        # so using file output instead
        tmp_blast_output_file = os.path.join("blastp.results" + str(uuid.uuid4()))

        blastp_cmd = f'blastp -outfmt 6 -subject {target_seq_file} -num_threads {threads} ' \
                     f'-query {query_seq_file} > {tmp_blast_output_file}'
        os.system(blastp_cmd)

        top_scores = [0] * number_vals_to_report
        top_records = [[]] * number_vals_to_report

        with open(tmp_blast_output_file) as bl:
            for line in bl:
                cols = blast_record(*line.strip().split())
                bl_score = float(cols.bitscore)
                if bl_score > noise_level and bl_score > min(top_scores):
                    min_i = top_scores.index(min(top_scores))
                    top_scores[min_i] = bl_score
                    top_records[min_i] = cols._asdict()

        return sorted(top_records, reverse=True, key=lambda r: r['bitscore'])

    def find_genes_from_similar_reactions(self, params):
        self._validate_params(params, {'workspace_name', 'reaction_set', 'query_genome_ref', },
                              {'number_of_hits_to_report', 'smarts_set', 'blast_score_floor',
                               'structural_similarity_floor', 'difference_similarity_floor'})

        feature_seq_path = self.gfu.genome_proteins_to_fasta(
            {'genome_ref': params['query_genome_ref'],
             'include_functions': True,
             'include_aliases': False})['file_path']

        # TODO: fix UI so I don't have to use this hack
        if isinstance(params['reaction_set'], str):
            params['reaction_set'] = [params['reaction_set']]

        output = {'gene_hits': [self.find_genes_for_rxn(rxn, feature_seq_path, params)
                                for rxn in params['reaction_set']]}
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
                                        params.get('blast_score_floor', 50),
                                        params.get('number_of_hits_to_report', 5))
                                        #params.get('number_vals_to_report', 5))

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
        tables = []
        for rxn, hits in zip(reactions, gene_hits):
            tables.append(f'<h3 style="text-align: center">{rxn}</h3>')
            if not hits:
                tables.append("<h5>No hits found!</h5>")
                continue
            tables.append('<table class="table table-bordered table-striped">')
            header = "</td><td>".join(hits[0].keys())
            tables.append(f'\t<thead><tr><td>{header}</td></tr></thead>')
            tables.append('\t<tbody>')
            for row in hits:
                line = "</td><td>".join(row.values())
                tables.append(f'\t\t<tr><td>{line}</td></tr>')
            tables.append('\t</tbody>')
            tables.append('</table>\n')

        # Fill in template HTML
        with open(os.path.join(os.path.dirname(__file__), 'find_genes_for_rxn_template.html')
                  ) as report_template_file:
            report_template = report_template_file.read() \
                .replace('*TABLES*', "\n".join(tables))

        with open(result_file_path, 'w') as result_file:
            result_file.write(report_template)

        html_report.append({'path': output_directory,
                            'name': os.path.basename(result_file_path),
                            'description': 'HTML report for reaction to gene candidates app'})

        return html_report
