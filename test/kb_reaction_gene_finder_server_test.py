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
        self.assertLessEqual( r, 0.04 )  # four percent tolerance

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

    #def test_run_kb_reaction_gene_finder_3(self):
    #
    #    inp = {'workspace_name': self.wsName,
    #           'reaction_set': ['rxn39253'],     # not found in modelseed, no similars
    #           'structural_similarity_floor': 1.0,
    #           'difference_similarity_floor': 1.0,
    #           'query_genome_ref': 'ReferenceDataManager/GCF_000009625.1',
    #           'number_of_hits_to_report': 5
    #           }
    #
    #    ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
    #             
    #    pprint(ret)
    #    self.validateRetStruct( inp, ret )
    #    self.validateValues( inp, ret, [[]] )


    # this case is problematic - the differences between blast scores go up to 3.5% and
    # the ordering throws off the comparison beyond the top 3 hits.

    #def test_run_kb_reaction_gene_finder_4(self):
    #
    #    inp = {'workspace_name': self.wsName,
    #           'reaction_set': ['rxn14379'],
    #           'structural_similarity_floor': 0.9,
    #           'difference_similarity_floor': 0.9,
    #           'query_genome_ref': 'ReferenceDataManager/GCF_001678945.1',
    #           'number_of_hits_to_report': 4
    #           }
    #
    #    ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
    #             
    #    pprint(ret)
    #    self.validateRetStruct( inp, ret )
    #    self.validateValues( inp, ret,
    #        [[['JL2886_RS04030', 'SL003B_RS05585', '59.531', '341', '137', '1', '13', '353', '13', '352', '5.98e-138', '387'],
    #          ['JL2886_RS04030', 'RL_RS33440', '60.882', '340', '132', '1', '14', '353', '14', '352', '1.38e-136', '384'],
    #          ['JL2886_RS04030', 'SADFL11_RS12560', '57.507', '353', '149', '1', '1', '353', '1', '352', '9.89e-136', '382'],
    #          #['JL2886_RS04030', 'H009_RS06465', '59.118', '340', '138', '1', '13', '352', '12', '350', '4.20e-130', '367'],
    #          #['JL2886_RS04030', 'VEIS_RS22320', '55.000', '340', '153', '0', '13', '352', '12', '351', '4.42e-122', '347'],
    #          #['JL2886_RS04030', 'YUQ_RS0102485', '53.529', '340', '158', '0', '13', '352', '13', '352', '1.17e-118', '339'],
    #          #['JL2886_RS04030', 'REUT_RS08005', '50.872', '344', '164', '2', '13', '352', '12', '354', '1.04e-116', '333'],
    #          ['JL2886_RS04030', 'P350_RS12650', '53.529', '340', '158', '0', '13', '352', '11', '350', '4.22e-116', '332']
    #          #['JL2886_RS04030', 'REUT_RS23880', '49.709', '344', '173', '0', '9', '352', '8', '351', '1.71e-109', '315'],
    #          #['JL2886_RS04030', 'PLU_RS00725', '46.784', '342', '180', '1', '13', '352', '12', '353', '2.50e-109', '315']
    #         ]] )

    def test_run_kb_reaction_gene_finder_5(self):
    
        inp = {'workspace_name': self.wsName,
               'reaction_set': ['rxn02432'],        # singlular occurrence in ModelSEED, no equivalents
               'structural_similarity_floor': 1.0,
               'difference_similarity_floor': 1.0,
               'query_genome_ref': 'ReferenceDataManager/GCF_000016765.1',  # occurs in this genome
               'number_of_hits_to_report': 10   # should not get this many
               }
    
        ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
                 
        pprint(ret)
        self.validateRetStruct( inp, ret )
        self.validateValues( inp, ret,
            [[
             ['SWIT_RS04250','SWIT_RS04250','100.000','793','0','0','1','793','1','793','0.0','1577'],
             ['SWIT_RS02535','SWIT_RS04250','50.066','753','357','6','15','753','23','770','0.0','630'],
             ['SWIT_RS23205','SWIT_RS04250','','','','','','','','','','154'],  # wasn't in offline run

             #['SWIT_RS23205','SWIT_RS04250','44.211','95','50','3','358','449','161','255','3.53e-12','56.2'],
             #['SWIT_RS25760', 'SWIT_RS04250','35.211', '142','88','3','608','747','646','785','3.33e-21','85.5'],
             #['SWIT_RS25760', 'SWIT_RS04250', '23.608','521', '350','11','23','505','25','535','1.16e-19', '80.5']
            ]] )

           
