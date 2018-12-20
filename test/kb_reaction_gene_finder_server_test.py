# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser
from pprint import pprint

from kb_reaction_gene_finder.kb_reaction_gene_finderImpl import kb_reaction_gene_finder
from kb_reaction_gene_finder.kb_reaction_gene_finderServer import MethodContext
from kb_reaction_gene_finder.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class kb_reaction_gene_finderTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_reaction_gene_finder'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_reaction_gene_finder',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_reaction_gene_finder(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.prepareTestData()

    @classmethod
    def prepareTestData(cls):
        """This function creates testing objects"""
        pass

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def validateRetStruct( self,  inp, ret ):
        print( "#length of ret is {0}".format( len(ret) ) )
        print( "#length of hits: {0}". format( len(ret[0].get("gene_hits",[] ))))
        print( "#length of reaction set: {0}". format( len(inp.get( 'reaction_set' )) ) )
        self.assertEqual( len(inp.get( 'reaction_set' )), len(ret[0].get("gene_hits",[]) ) )
        for subl in ret[0].get( "gene_hits", [] ):
            print( "#length of subl is {0}".format( len( subl ) ) )
            self.assertLessEqual( len( subl ), inp.get( 'number_of_hits_to_report', 5 ) )
        print( "#return structure validated" )

    def validateRow( self, vals, e_vals ):
        print( "#validate rows {0} vs {1}".format( vals, e_vals ) )
        print( "compare {0} vs {1}".format( vals.get( 'qseqid'), e_vals[1] ) )
        self.assertEqual( vals.get( 'qseqid' ), e_vals[1] )
        print( "compare {0} vs {1}".format( vals.get( 'bitscore'), e_vals[11] ) )
        r = abs( float( vals.get( 'bitscore' ) ) - float( e_vals[11] ) ) / float( e_vals[11] )
        self.assertLessEqual( r, 0.02 )  # two percent tolerance

    def validateValues( self, inp, ret, expected_vals ):
        n = len(ret)
        #n = len(ret[0].get("gene_hits",[]) )
        self.assertEqual( len( expected_vals ), n )
        for i in range(0,n):
            print( "#checking list {0}:".format( i ) )
            vals = ret[0].get("gene_hits")[i]
            e_vals = expected_vals[i]
            nrows = len( e_vals )
            print( "#len vals {0} vs e_vals {1}".format( len( vals ), len( e_vals ) ) )
            self.assertEqual( nrows, len( vals ) )
            for j in range(0,nrows):
                print( "#checking row {0}:".format( j ) )
                self.validateRow( vals[j], e_vals[j] )
        print( "#values validated" )
                

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    #def test_run_kb_reaction_gene_finder_ok(self):
    #    # call your implementation
    #    inp = {'workspace_name': self.wsName,
    #           'reaction_set': 'rxn00010',
    #           'query_genome_ref': 'ReferenceDataManager/GCF_002163935.1',
    #           'number_of_hits_to_report': 10
    #           }
    #    ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
    #    pprint(ret)
    #    self.validateRetStruct( inp, ret )

    #def test_run_kb_reaction_gene_finder_2(self):
    #
    #    inp = {'workspace_name': self.wsName,
    #           'reaction_set': ['rxn00008'],     # this reaction is not in modelseeds, but similars are
    #           'structural_similarity_floor': 0.0,
    #           'difference_similarity_floor': 1.0,
    #           'query_genome_ref': 'ReferenceDataManager/GCF_000009625.1',
    #           'number_of_hits_to_report': 5
    #           }
    #
    #    ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
    #             
    #    pprint(ret)
    #    self.validateRetStruct( inp, ret )
    #    self.validateValues( inp, ret, 
    #      [[['MAFF_RS23790', 'MAFF_RS23790', '100.000', '361', '0', '0', '1', '361', '1', '361', '0.0', '743'],
    #        ['MAFF_RS23790', 'DSHI_RS13475', '71.605', '324', '92', '0', '37', '360', '27', '350', '0.0', '508'],
    #        ['MAFF_RS23790', 'ROS217_RS07690', '69.940', '336', '101', '0', '23', '358', '4', '339', '0.0', '504'],
    #        ['MAFF_RS23790', 'COC_RS0105080', '72.050', '322', '90', '0', '37', '358', '23', '344', '2.74e-180', '498'],
    #        ['MAFF_RS23790', 'RSPH17029_RS05395', '71.739', '322', '91', '0', '37', '358', '23', '344', '4.62e-179', '495']]]
    #       )

    def test_run_kb_reaction_gene_finder_3(self):
    
        inp = {'workspace_name': self.wsName,
               'reaction_set': ['rxn39253'],     # not found in modelseed, no similars
               'structural_similarity_floor': 1.0,
               'difference_similarity_floor': 1.0,
               'query_genome_ref': 'ReferenceDataManager/GCF_000009625.1',
               'number_of_hits_to_report': 5
               }
    
        ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
                 
        pprint(ret)
        self.validateRetStruct( inp, ret )
        self.validateValues( inp, ret, [[]] )

           
