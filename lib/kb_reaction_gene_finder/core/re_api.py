import json
import logging
from pprint import pformat

import requests


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
        if not rid.startswith('rxn_reaction/'):
            rid = 'rxn_reaction/' + rid
        body = json.dumps({'rid': rid, 'sf_sim': sf_sim, 'df_sim': df_sim,
                           'exclude_self': exclude_self})
        ret = self._call_re(params={'view': "list_genes_for_similar_reactions"}, data=body)
        if "error" in ret:
            raise RuntimeError(f"{ret['error']}: {ret.get('arango_message', '')}")
        logging.info(f"Found {ret['results'][0]['count']} related sequences")
        return ret['results'][0]['sequences']

