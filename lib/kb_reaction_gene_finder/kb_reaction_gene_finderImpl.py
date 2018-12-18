# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import logging
import os
from pprint import pformat

from Bio import SeqIO

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class reaction_gene_finder:
    '''
    Module Name:
    reaction_gene_finder

    Module Description:
    note on terms:  "query" means the unannotaed genes in the users unannotated genome
"subject" means a matching gene from the modelseed set or a
matching uniref sequence mapped from the modelseed set
 I
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/janakagithub/kb_rxn_gene_finder.git"
    GIT_COMMIT_HASH = "b3b49759b73bb1f5b251c8e90b661179d7c326cc"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        
        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def find_genes_for_exact_rxn_matches(self, ctx, params):
        """
        :param params: instance of type "find_genes_exact_matches" ->
           structure: parameter "workspace" of String, parameter
           "reaction_set" of list of String, parameter "smarts_set" of list
           of String, parameter "query_genome_ref" of String, parameter
           "number_of_hits_to_report" of Long
        :returns: instance of type "find_genes_exact_matches_results" ->
           structure: parameter "gene_hits_info" of list of type "gene_hits"
           -> structure: parameter "reaction_id" of String, parameter
           "smarts_id" of String, parameter "structural_similarity_score" of
           Double, parameter "difference_similarity_score" of Double,
           parameter "structural_distance" of Double, parameter
           "difference_distance" of Double, parameter "top_gene_hits" of
           mapping from String to list of String, parameter "reactions" of
           list of String, parameter "query_genome_ref" of String, parameter
           "query_genome_name" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN find_genes_for_exact_rxn_matches
        #END find_genes_for_exact_rxn_matches

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method find_genes_for_exact_rxn_matches return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
