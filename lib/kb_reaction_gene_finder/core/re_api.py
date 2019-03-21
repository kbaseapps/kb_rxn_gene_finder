import json
import logging
from pprint import pformat

import requests


class RE_API:
    def __init__(self, re_url, token):
        self.re_url = re_url
        self.token = token

    def _call_re(self, endpoint="/api/v1/query_results/", params=None, data=None):
        header = {"Authorization": self.token}
        logging.info(f"Calling RE_API with query data: {pformat(data)}")
        ret = requests.post(self.re_url+endpoint, data, params=params, headers=header)
        return ret.json()

    def get_related_sequences_adhoc(self, rid, sf_sim=1, df_sim=1, exclude_self=0):
        query = """
        WITH rxn_reaction
        LET ws_ids = @ws_ids
        LET start = @exclude_self ? 1 : 0
        LET rxns = (
            FOR v, e IN start..1
                ANY @rid rxn_similar_to_reaction
                OPTIONS {uniqueVertices: "global", bfs: true}
                FILTER !e || e.sf_similarity >= @sf_sim
                FILTER !e || e.df_similarity >= @df_sim
                RETURN {id: v._id, key: v._key, name: v.name, definition: v.definition, "structural similarity": e.sf_similarity, "difference similarity": e.df_similarity}
        )
        LET rxn_ids = rxns[*].id
        
        LET rxn_gene_links = (
            FOR e in rxn_reaction_within_complex
                FILTER e._from in rxn_ids
                LET linked_gene_ids = FLATTEN(
                    FOR c in rxn_gene_complex
                       FILTER c._id == e._to
                       RETURN c.genes
                )
                COLLECT rxn_id = e._from INTO groups KEEP linked_gene_ids
                RETURN {rxn_id: rxn_id, linked_gene_ids: UNIQUE(FLATTEN(groups[*].linked_gene_ids))}
        )
        
        LET gene_ids = UNIQUE(FLATTEN(rxn_gene_links[*].linked_gene_ids))
        
        LET genes = (
            FOR g in ncbi_gene
               FILTER g._key IN gene_ids
               RETURN {key: g._key, product: g.product, function: CONCAT_SEPARATOR(', ', g.functions), sequence: g.protein_translation}
        )
        LET missing = MINUS(gene_ids, genes[*].key)
        
        RETURN {rxns: rxns, rxn_gene_links: rxn_gene_links, genes: genes, missing_genes: missing}
        """
        if not rid.startswith('rxn_reaction/'):
            rid = 'rxn_reaction/' + rid
        body = json.dumps({'rid': rid, 'sf_sim': sf_sim, 'df_sim': df_sim,
                           'exclude_self': exclude_self, "query": query})
        ret = self._call_re(data=body)
        if "error" in ret:
            raise RuntimeError(f"{ret['error']}: {ret.get('arango_message', '')}")
        logging.info(f"Found {len(ret['results'][0]['genes'])} related sequences")
        return ret['results'][0]

    def get_related_sequences(self, rid, sf_sim=1, df_sim=1, exclude_self=0, retries=2):
        if not rid.startswith('rxn_reaction/'):
            rid = 'rxn_reaction/' + rid
        body = json.dumps({'rid': rid, 'sf_sim': sf_sim, 'df_sim': df_sim,
                           'exclude_self': exclude_self})
        ret = self._call_re(params={'view': "list_genes_for_similar_reactions"}, data=body)
        if "error" in ret:
            if retries:
                logging.warning("Arango Query failed. Retrying")
                return self.get_related_sequences(rid, sf_sim, df_sim, exclude_self, retries-1)
            raise RuntimeError(f"{ret['error']}: {ret.get('arango_message', '')}")
        logging.info(f"Found {len(ret['results'][0]['genes'])} related sequences")
        return ret['results'][0]

