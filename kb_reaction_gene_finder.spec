/*
A KBase module: Gene predictions based on finger print reactions
This module predict potential gene candidates based on exact reactions or reaction fingerprints
*/

/* note on terms:  "query" means the unannotaed genes in the users unannotated genome
                   "subject" means a matching gene from the modelseed set or a
                   matching uniref sequence mapped from the modelseed set
                    I
*/
module kb_reaction_gene_finder {

    typedef structure {
        string workspace_name;
        list <string> reaction_set;
    	list <string> smarts_set;
        string query_genome_ref;
        float structural_similarity_floor;
        float difference_similarity_floor;
        float blast_score_floor;
        int number_of_hits_to_report;
    } findGenesParams;

	/*
     for non-exact matches, we would probably add at least two extra parameters,
     which would be cutoff thresholds for the two similarity scores
	*/

    typedef structure{
        string query_gene_id;          /* this is users "unannotated" gene */
        string subject_gene_id;        /* this is the matching modelseed/KEGG gene id */
        string subject_genome;         /* this is the matching modelseed/KEGG genome id */
        string subject_genome_ref;	   /* this is the matching modelseed/KEGG genome ref */
        string gene_db;                /* field for "Modelseed" or "KEGG" */

		/*
		I'm adding in an option for using UniRef50 or not
        string Uniref cluster_id blank if not using uniref */

        string uniref50_cluster_id;

        float  blast_score;            /* these four come from the blast report*/
        float  percent_identity;
        float  match_length;
        float  coverage;
        float  frequency;    /* Frequency as to how common this gene occurs in genomes */
    } top_gene_hits;


    typedef structure{
       string reaction_id;
       string smarts_id;
       float structural_similarity_score;
       float difference_similarity_score;
       mapping <string, list <string>> top_gene_hits;
    } GeneHits;


    typedef structure {

        list <GeneHits> gene_hits;

    } findGenesResults;


    funcdef find_genes_from_similar_reactions(findGenesParams params) returns (findGenesResults output)
        authentication required;
};
