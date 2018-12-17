/*
A KBase module: Gene predictions based on finger print reactions
This module predict potential gene candidates based on exact reactions or reaction fingerprints
*/

/* note on terms:  "query" means the unannotaed genes in the users unannotated genome
                   "subject" means a matching gene from the modelseed set or a
                   matching uniref sequence mapped from the modelseed set
                    I
*/
module reaction_gene_finder {

    typedef structure {
        string workspace;
        list <string> reaction_set;
    	list <string> smarts_set;
        string query_genome_ref;
        int number_of_hits_to_report;   /* default 5?  10? */
    } find_genes_exact_matches;

	/*
     for non-exact matches, we would probably add at least two extra parameters,
     which would be cutoff thresholds for the two similarity scores
	*/
    typedef structure{
        string query_gene_id;          /* this is users "unannotated" gene */
        string subject_gene_id;        /* this is the matching modelseed gene id */
        string subject_genome;         /* this is the matching modelseed genome id */
        /*
         maybe we should include a flag or fieldfor "Modelseed" or "KEGG"?
        /*
        string uniref50_cluster_id;    /* blank if not using uniref */

        float  blast_score;            /* these four come from the blast report*/
        float  percent_identity;
        float  match_length;
        float  coverage;
        float  frequency;    /* ? not sure what this is, left it in*/
    } top_gene_hits;


    typedef structure{
        /* I'm adding in an option for using UniRef50 or not
        string Uniref cluster_id
        */
        string genome_id;
        mapping <string, list <string>> gene_ids_scores;
    } gene_ids_by_genome;


    typedef structure{
       string reaction_id;
       string smarts_id;
       float structural_similarity_score;
       float difference_similarity_score;
       float structural_distance;     /* (1.0 - structural_similarity_score) */
       float difference_distance;     /* (1.0 - difference_similarity_score)*/
       /* note: for this exact case, both similarities will be 1 and
        both distances will be 0
       */
       mapping <string, list <string>> top_gene_hits;
    } gene_hits;


    typedef structure {

        list <gene_hits> gene_hits_info;

        /*Any metadata from the query
        how many genomes
        Uniref info
        */

        list <string> reactions;
        string query_genome_ref;
        string query_genome_name;

    } find_genes_exact_matches_results;


    funcdef find_genes_from_exact_matches(find_genes_exact_matches params) returns (find_genes_exact_matches_results output)
        authentication required;
};
