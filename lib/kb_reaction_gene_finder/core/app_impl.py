import logging
import os
import uuid
from collections import OrderedDict, Counter

from kb_reaction_gene_finder.core.re_api import RE_API
from installed_clients.FeatureSetUtilsClient import FeatureSetUtils
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.KBaseReportClient import KBaseReport


def _make_table_html(title, items, keys=None, empty_message="No items found"):
    """Takes a list of dicts and makes a HTML table"""
    table_lines = []
    table_lines.append(f'<h4 style="text-align: center">{title}</h4>')
    if items:
        if not keys:
            keys = items[0].keys()
        table_lines.append('<table class="table table-bordered table-striped">')
        header = "</td><td>".join(x.capitalize() for x in keys)
        table_lines.append(f'\t<thead><tr><td>{header}</td></tr></thead>')
        table_lines.append('\t<tbody>')
        for row in items:
            line = "</td><td>".join((str(row[k]) if row.get(k) is not None else "" for k in keys))
            table_lines.append(f'\t\t<tr><td>{line}</td></tr>')
        table_lines.append('\t</tbody>')
        table_lines.append('</table>\n')
    else:
        table_lines.append(f"<p>{empty_message}<p>")
    return "\n".join(table_lines)


def _make_rxn_html(arango_results, gene_hits):
    rxn_tbl = _make_table_html("Similar Reactions", arango_results.get('rxns'),
                               ('key', 'name', 'definition', 'structural similarity',
                                'difference similarity'),
                               "No similar reactions were found. You may need to relax "
                               "the Structural or Difference Similarity Floors.")

    gene_tbl = _make_table_html("Related Genes", arango_results.get('genes'),
                                ('key', 'product', 'function'),
                                "No genes related to the similar reactions above were found in our"
                                " database. You may need to relax the Structural or Difference "
                                "Similarity Floors.")
    hits_tbl = _make_table_html("Genome Hits", gene_hits,
                                empty_message="No related genes were close matches to your genome."
                                              " If you have 'Related Genes' you may need to change"
                                              " the Blast Score Floor.")
    return "\n".join([rxn_tbl, gene_tbl, hits_tbl])


class AppImpl:
    def __init__(self, config, ctx):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.re_api = RE_API(config['re-api'], ctx['token'])
        self.fsu = FeatureSetUtils(self.callback_url)
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

        if not params.get('reaction_set'):
            params['reaction_set'] = []

        if params.get('bulk_reaction_ids'):
            params['reaction_set'] += params['bulk_reaction_ids'].split('\n')

        if not params['reaction_set']:
            raise ValueError("No reactions to analyze")

        return params['reaction_set']

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

    def _make_feature_set(self, workspace, genome, fs_name_prefix, reaction_id, top_genes):
        params = {'genome': genome,
                  'feature_ids': top_genes,
                  'workspace_name': workspace,
                  'description': f'A set of the top gene candidates for {reaction_id} calculated '
                                 f'by the "Find Candidate Genes for a Reaction" app',
                  'output_feature_set': f'{fs_name_prefix}_{reaction_id}',
                 }
        return self.fsu.build_feature_set(params)['feature_set_ref']

    def _find_best_homologs(self, query_seq_file, target_seq_file, rxn_id,
                            noise_level=50, number_vals_to_report=5, threads=1):
        """Blast the query_seq_file against the target_seq_file and return the best hits"""
        logging.info("running blastp for {0} vs {1}".format(query_seq_file, target_seq_file))
        out_cols = ['sseqid', 'qseqid', 'bitscore', 'pident', 'length', 'mismatch', 'evalue',]
        col_names = ['Genome Gene', 'Closest Database Gene', 'Bit Score', 'Percent Identity',
                     'Match Length', 'Mismatches', 'E Value']

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

        top_genes = [gene[0] for gene in top_bitscore.most_common(number_vals_to_report)]
        top_records = [{"Reaction ID": rxn_id, **gene_hits[gene],
                        "Total Gene Hits": str(gene_hit_count[gene])}
                       for gene in top_genes]

        return top_records, top_genes

    def find_genes_from_similar_reactions(self, params):
        reaction_ids = self._validate_params(
            params, {'workspace_name', 'query_genome_ref', },
            {'number_of_hits_to_report', 'smarts_set', 'blast_score_floor',
            'structural_similarity_floor', 'difference_similarity_floor',
            'reaction_set', 'bulk_reaction_ids'})

        feature_seq_path = self.gfu.genome_proteins_to_fasta(
            {'genome_ref': params['query_genome_ref'],
             'include_functions': True,
             'include_aliases': False})['file_path']

        output = {'gene_hits': [], 'feature_set_refs': []}
        html_tables = []
        for rxn in reaction_ids:
            hits, genes, html = self.find_genes_for_rxn(rxn, feature_seq_path, params)
            if genes:
                output['feature_set_refs'].append(
                    self._make_feature_set(params['workspace_name'],
                                           params['query_genome_ref'],
                                           params.get('feature_set_prefix', 'gene_candidates'),
                                           rxn,
                                           genes))
            output['gene_hits'].extend(hits)
            print(html)
            html_tables.append(html)
        output.update(self._build_report(reaction_ids,
                                         html_tables,
                                         output['feature_set_refs'],
                                         params['workspace_name'],
                                         ))
        return output

    def find_genes_for_rxn(self, reaction, genome_feature_path, params):
        """Finds genes for a particular reaction using RE and BLAST"""
        arango_results = self.re_api.get_related_sequences_adhoc(
            reaction,
            params.get('structural_similarity_floor', 1),
            params.get('difference_similarity_floor', 1))
        if not arango_results.get('genes'):
            return [], [], _make_rxn_html(arango_results, [])

        search_fasta = f'{self.scratch}/rxn_search_{uuid.uuid4()}.fasta'
        self._make_fasta(arango_results['genes'], search_fasta)

        hits, genes = self._find_best_homologs(search_fasta,
                                        genome_feature_path,
                                        reaction,
                                        params.get('blast_score_floor', 50),
                                        params.get('number_of_hits_to_report', 5))
        html = _make_rxn_html(arango_results, hits)
        return hits, genes, html

    def _build_report(self, reaction_ids, html_tables, feature_sets, workspace_name):
        """
        _generate_report: generate summary report for upload
        """
        output_html_files = self._generate_report_html(reaction_ids, html_tables)

        report_params = {
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'objects_created': [{'ref': fs_ref, 'description': 'A set of the top gene candidates'}
                                for fs_ref in feature_sets],
            'workspace_name': workspace_name,
            'report_object_name': 'find_genes_for_rxn_' + str(uuid.uuid4())}

        output = self.kbr.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_report_html(self, reactions, html_tables):
        """
            _generate_report: generates the HTML for the upload report
        """
        assert len(reactions) == len(html_tables)

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'find_genes_for_rxn.html')

        # first line is special (because it should be open)
        first_id = reactions.pop(0)
        lines = ['<div class="tab">',
                 f'\t<button id="defaultOpen" class="tablinks" '
                 f'onclick="openTab(event, \'{first_id}\')">{first_id}</button>']

        lines += [f'\t<button class="tablinks" onclick="openTab(event, \'{rid}\')">{rid}</button>'
                  for rid in reactions]
        lines.append('</div>')
        lines += [f'<div id={rid} class="tabcontent">\n{html}\n</div>'
                  for rid, html in zip([first_id]+reactions, html_tables)]

        # Fill in template HTML
        with open(os.path.join(os.path.dirname(__file__), 'find_genes_for_rxn_template.html')
                  ) as report_template_file:
            report_template = report_template_file.read() \
                .replace('*TABLES*', "\n".join(lines))

        with open(result_file_path, 'w') as result_file:
            result_file.write(report_template)

        html_report = [{'path': output_directory,
                        'name': os.path.basename(result_file_path),
                        'description': 'HTML report for reaction to gene candidates app'}]

        return html_report
