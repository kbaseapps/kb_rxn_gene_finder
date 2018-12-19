import logging
import os
import time

from kb_reaction_gene_finder.core import app_skeleton
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

    def find_genes_from_exact_matches(self, params):
        self._validate_params(params, {'workspace_name', 'reaction_set', 'query_genome_ref', },
                              {'number_of_hits_to_report', 'smarts_set'})

        feature_seq_path = self.gfu.genome_proteins_to_fasta(
            {'genome_ref': params['query_genome_ref'],
             'include_functions': False,
             'include_aliases': False})['file_path']
        # TODO: add in a report
        return {'gene_hits': [self.find_genes_for_rxn(rxn, feature_seq_path, params)
                              for rxn in params['reaction_set']]}

    def find_genes_for_rxn(self, reaction, genome_feature_path, params):
        search_seqs = self.re_api.get_related_sequences(reaction)
        if not search_seqs:
            return "None found"
        search_fasta = f'rxn_search_{time.time()}.fasta'
        self._make_fasta(search_seqs, search_fasta)
        return app_skeleton.find_best_homologs(search_fasta, genome_feature_path)
