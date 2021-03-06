# -*- coding: utf-8 -*-
import re
import os
import time
import unittest
from configparser import ConfigParser
from collections import Counter
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

    def isGoodWSRef(self, s):
        if not isinstance(s, str):
            return (False)
        return (bool(re.match("^\d+\/\d+\/\d+$", s)))

    def checkMembers(self, d, memlist):
        for m in memlist:
            self.assertIn(m, d, "{0} not found in {1}".format(m, d))

    def validateRetStruct(self, inp, ret):
        print("validateRetStruct, inp is")
        pprint(inp)
        print("ret is ")
        pprint(ret)

        # Make sure it is a list of structures with the correct members

        self.assertIsInstance(ret, list, "return must be a list")
        for r in ret:
            self.checkMembers(r, ['gene_hits', 'feature_set_refs', 'report_name', 'report_ref'])
            self.assertIsInstance(r.get('report_name'), str,
                                  "report_name is not string {0}".format(r))
            self.assertTrue(self.isGoodWSRef(r.get('report_ref')),
                            f"report_ref is not a good workspace object ref: {r}")
            fl = r.get('feature_set_refs')
            self.assertIsInstance(fl, list,
                                  f"feature_set_refs {fl} is not list ({fl})")
            for f in fl:
                self.assertTrue(self.isGoodWSRef(f),
                                f"{f} is not a good ws object ref ({fl})")
            gl = r.get('gene_hits')
            self.assertIsInstance(gl, list,
                                  f"gene_hits {gl} is not list ({gl})")
            for g in gl:
                self.assertIsInstance(g, dict, f"gene hit {g} is not dict ({gl})")
                self.checkMembers(g, ['Bit Score', 'Closest Database Gene', 'E Value',
                                      'Genome Gene', 'Match Length', 'Mismatches',
                                      'Percent Identity', 'Total Gene Hits'])

        print("#return structure validated")

    def validateRow(self, vals, e_vals):
        self.assertEqual(vals.get('Database Gene'), e_vals[1])
        r = abs(float(vals.get('Bit Score')) - float(e_vals[-1])) / float(e_vals[-1])
        self.assertLessEqual(r, 0.04)  # four percent tolerance

    def validateValues(self, inp, ret, expected_vals):
        n = len(ret)
        self.assertEqual(len(expected_vals), n)
        for i in range(0, n):
            vals = ret[0].get("gene_hits")[i]
            e_vals = expected_vals[i]
            nrows = len(e_vals)
            self.assertEqual(nrows, len(vals))
            for j in range(0, nrows):
                self.validateRow(vals[j], e_vals[j])

    def test_find_genes_from_similar_reactions_bad_input(self):
        with self.assertRaisesRegex(ValueError, "No reactions to analyze"):
            inp = {'workspace_name': self.wsName,
                   'query_genome_ref': 'ReferenceDataManager/GCF_002163935.1',
                   'number_of_hits_to_report': 10
                   }
            ret = self.serviceImpl.find_genes_from_similar_reactions(self.ctx, inp)
        with self.assertRaisesRegex(ValueError, "workspace_name not in supplied parameters"):
            inp = {'reaction_set': ['rxn00008'],
                   'query_genome_ref': 'ReferenceDataManager/GCF_002163935.1',
                   'number_of_hits_to_report': 10
                   }
            ret = self.serviceImpl.find_genes_from_similar_reactions(self.ctx, inp)
        with self.assertRaisesRegex(ValueError, "query_genome_ref not in supplied parameters"):
            inp = {'reaction_set': ['rxn00008'],
                   'workspace_name': self.wsName,
                   'number_of_hits_to_report': 10
                   }
            ret = self.serviceImpl.find_genes_from_similar_reactions(self.ctx, inp)

    def test_find_genes_from_similar_reactions_no_hits(self):
        inp = {'workspace_name': self.wsName,
               'reaction_set': ['rxn01965'],
               'query_genome_ref': 'ReferenceDataManager/GCF_002163935.1',
               'number_of_hits_to_report': 10
               }
        ret = self.serviceImpl.find_genes_from_similar_reactions(self.ctx, inp)

    def test_find_genes_from_similar_reactions_string(self):
        inp = {'workspace_name': self.wsName,
               'bulk_reaction_ids': 'rxn00010\nrxn14379',
               'query_genome_ref': 'ReferenceDataManager/GCF_002163935.1',
               'structural_similarity_floor': 0.7,
               'difference_similarity_floor': 0.7,
               'number_of_hits_to_report': 10
               }
        ret = self.serviceImpl.find_genes_from_similar_reactions(self.ctx, inp)
        self.validateRetStruct(inp, ret)

    def test_find_genes_from_similar_reactions_2rxn_by10(self):
        inp = {'workspace_name': self.wsName,
               'bulk_reaction_ids': 'rxn00371\nrxn00083\nrxn04632',
               'query_genome_ref': 'ReferenceDataManager/GCF_002163935.1',
               'number_of_hits_to_report': 10
               }
        ret = self.serviceImpl.find_genes_from_similar_reactions(self.ctx, inp)
        self.validateRetStruct(inp, ret)

    # return value checks

    """
    def test_find_genes_from_similar_reactions_2(self):

       inp = {'workspace_name': self.wsName,
              'reaction_set': ['rxn00008'],   # this reaction is not in modelseeds, but similars are
              'structural_similarity_floor': 0.0,
              'difference_similarity_floor': 1.0,
              'query_genome_ref': 'ReferenceDataManager/GCF_000009625.1',
              'number_of_hits_to_report': 5
              }

       ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
       self.validateRetStruct( inp, ret )
       self.validateValues( inp, ret, 
         [[['MAFF_RS23790', 'MAFF_RS23790', '100.000', '361', '0', '0.0', '743'],
           ['MAFF_RS23790', 'DSHI_RS13475', '71.605', '324', '92', '0.0', '508'],
           ['MAFF_RS23790', 'ROS217_RS07690', '69.940', '336', '101', '0.0', '504'],
           ['MAFF_RS23790', 'COC_RS0105080', '72.050', '322', '90', '2.74e-180', '498'],
           ['MAFF_RS23790', 'RSPH17029_RS05395', '71.739', '322', '91', '4.62e-179', '495']]]
          )

    def test_find_genes_from_similar_reactions_3(self):

       inp = {'workspace_name': self.wsName,
              'reaction_set': ['rxn39253'],     # not found in modelseed, no similars
              'structural_similarity_floor': 1.0,
              'difference_similarity_floor': 1.0,
              'query_genome_ref': 'ReferenceDataManager/GCF_000009625.1',
              'number_of_hits_to_report': 5
              }

       ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
       self.validateRetStruct( inp, ret )
       self.validateValues( inp, ret, [[]] )

    @unittest.skip("this case is problematic - the differences between blast scores go up to 3.5% "
                   "and the ordering throws off the comparison beyond the top 3 hits.")
    def test_find_genes_from_similar_reactions_4(self):

       inp = {'workspace_name': self.wsName,
              'reaction_set': ['rxn14379'],
              'structural_similarity_floor': 0.9,
              'difference_similarity_floor': 0.9,
              'query_genome_ref': 'ReferenceDataManager/GCF_001678945.1',
              'number_of_hits_to_report': 4
              }

       ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )

       pprint(ret)
       self.validateRetStruct( inp, ret )
       self.validateValues( inp, ret,
           [[['JL2886_RS04030', 'SL003B_RS05585', '59.531', '341', '137', '5.98e-138', '387'],
             ['JL2886_RS04030', 'RL_RS33440', '60.882', '340', '132', '1.38e-136', '384'],
             ['JL2886_RS04030', 'SADFL11_RS12560', '57.507', '353', '9.89e-136', '382'],
             ['JL2886_RS04030', 'P350_RS12650', '53.529', '340', '158', '4.22e-116', '332']
            ]] )

    def test_find_genes_from_similar_reactions_5(self):
    
        inp = {'workspace_name': self.wsName,
               'reaction_set': ['rxn02432'],    # singlular occurrence in ModelSEED, no equivalents
               'structural_similarity_floor': 1.0,
               'difference_similarity_floor': 1.0,
               'query_genome_ref': 'ReferenceDataManager/GCF_000016765.1',  # occurs in this genome
               'number_of_hits_to_report': 10   # should not get this many
               }
    
        ret = self.serviceImpl.find_genes_from_similar_reactions( self.ctx, inp )
        self.validateRetStruct( inp, ret )
        self.validateValues( inp, ret,
            [[
             ['SWIT_RS04250','SWIT_RS04250','100.000','793','0', '0.0','1577'],
             ['SWIT_RS02535','SWIT_RS04250','50.066','753','357', '0.0','630'],
             ['SWIT_RS23205','SWIT_RS04250','','','', '','154'],  # wasn't in offline run

             #['SWIT_RS23205','SWIT_RS04250','44.211','95','50','3','358','449','161','255','3.53e-12','56.2'],
             #['SWIT_RS25760', 'SWIT_RS04250','35.211', '142','88','3','608','747','646','785','3.33e-21','85.5'],
             #['SWIT_RS25760', 'SWIT_RS04250', '23.608','521', '350','11','23','505','25','535','1.16e-19', '80.5']
            ]] )"""
