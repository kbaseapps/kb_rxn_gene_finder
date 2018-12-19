# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import logging
import os

from kb_reaction_gene_finder.core.app_impl import AppImpl
#END_HEADER


class kb_reaction_gene_finder:
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
    GIT_URL = "https://github.com/kbaseapps/kb_rxn_gene_finder.git"
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
        self.config = config
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def find_genes_for_exact_rxn_matches(self, ctx, params):
        """
        :param params: instance of type "findGenesExactMatchesParams" ->
           structure: parameter "workspace_name" of String, parameter
           "reaction_set" of list of String, parameter "smarts_set" of list
           of String, parameter "query_genome_ref" of String, parameter
           "number_of_hits_to_report" of Long
        :returns: instance of type "findGenesExactMatchesResults" ->
           structure: parameter "gene_hits" of list of type "GeneHits" ->
           structure: parameter "reaction_id" of String, parameter
           "smarts_id" of String, parameter "structural_similarity_score" of
           Double, parameter "difference_similarity_score" of Double,
           parameter "top_gene_hits" of mapping from String to list of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN find_genes_from_exact_matches
        app_impl = AppImpl(self.config, ctx)
        output = app_impl.find_genes_from_exact_matches(params)
        #END find_genes_from_exact_matches

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
