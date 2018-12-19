import logging
import os
import uuid
from collections import namedtuple

from kb_reaction_gene_finder.core.re_api import RE_API
from installed_clients.GenomeFileUtilClient import GenomeFileUtil


class AppImpl:
    def __init__(self, config, ctx):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.re_api = RE_API(config['re-api'], ctx['token'])
        self.gfu = GenomeFileUtil(self.callback_url)

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
                            number_vals_to_report=5, threads=4):
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
             'include_functions': False,
             'include_aliases': False})['file_path']
        # TODO: add in a report
        return {'gene_hits': [self.find_genes_for_rxn(rxn, feature_seq_path, params)
                              for rxn in params['reaction_set']]}

    def find_genes_for_rxn(self, reaction, genome_feature_path, params):
        """Finds genes for a particular reaction using RE and BLAST"""
        search_seqs = self.re_api.get_related_sequences(
            reaction,
            params.get('structural_similarity_floor', 1),
            params.get('difference_similarity_floor', 1))
        if not search_seqs:
            return "None found"

        search_fasta = f'rxn_search_{uuid.uuid4()}.fasta'
        self._make_fasta(search_seqs, search_fasta)

        return self._find_best_homologs(search_fasta,
                                        genome_feature_path,
                                        params.get('number_vals_to_report', 5),
                                        params.get('blast_score_floor', 50))
