import json
import logging
from pprint import pformat

import requests

#
#  This is fapi - the fake api which is an attempt to isolate the layer of calls that
#  are anticipated to handled by the graphDB.   If this works, it should be possible to
#  rewrite these functions/procedures to use Arango and run the main app basis program
#  without alteration
#


class RE_API:
    def __init__(self, re_url, token):
        self.re_url = re_url
        self.token = token

    def _call_re(self, endpoint="/api/query_results/", params=None, data=None):
        header = {"Authorization": self.token}
        logging.info(f"Calling RE_API with query data: {pformat(data)}")
        ret = requests.post(self.re_url+endpoint, data, params=params, headers=header)
        return ret.json()

    def get_related_sequences(self, rid, sf_sim=1, df_sim=1, exclude_self=0):
        aql = '''
            let similar_rxn_ids = (
            for e in rxn_similar_to_reaction
               filter e.sf_similarity >= @sf_sim
               filter e.df_similarity >= @df_sim
               filter e._to == @rid || e._from == @rid
               return e._to == @rid ? e._from : e._to
            )
            let self = @exclude_self ? "no_self" : @rid 
            let similar_complex_ids = (
            for e in rxn_reaction_within_complex
               filter e._from in similar_rxn_ids || e._from == self
               return e._to
            )
            
            let genes = FLATTEN(
            for c in rxn_gene_complex
               filter c._id IN similar_complex_ids
               return c.genes
            )
            
            let sequences = (
            for g in ncbi_gene
               filter g._key IN genes
               return {key: g._key, seq: g.protein_translation}
            )
            
            return {count: COUNT(sequences), sequences: sequences}'''
        if not rid.startswith('rxn_reaction/'):
            rid = 'rxn_reaction/' + rid
        body = json.dumps({'rid': rid, 'sf_sim': sf_sim, 'df_sim': df_sim,
                           'exclude_self': exclude_self, 'query': aql})
        ret = self._call_re(data=body)
        logging.info(f"Found {ret['results'][0]['count']} related sequences")
        return ret['results'][0]['sequences']

