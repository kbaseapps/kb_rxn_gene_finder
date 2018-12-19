import os

from kb_reaction_gene_finder.core import fapi


#
# This assumes the following
#
#   1) reaction-reaction Jaccard similarities, from both structural and difference
#      fingerprints have been pre-calculated and are retrievable via a fake api call
#
#          fapi.get_similar_reactions()
#
#   2) modelSEED genome-gene
#
#   3) UniRef gene->UniRef clusters are available through a precalculated
#       map
#
#   4) a genome's proteome (complete set of all protein sequences) is available through a call
#
#         fapi.get_proteome( genome_id )
#

# this presumbably will be expanded to do a full-blown search of the reaction
# simlarity matrix, using thresholds and whatnot in "other_params"

def get_similar_reactions(reaction_id):
    # this similarity scan takes a couple of minutes to run
    # rxn_list = fapi.get_similar_reactions( reaction_id, other_params )
    # in the interests of speed, just bypassing it for now
    rxn_list = [reaction_id]
    return (rxn_list)


def get_model_genes_for_reaction(rxn):
    print("getting model genes for reaction {0}".format(rxn))
    return (fapi.get_modelseed_genomes_genes(rxn))


def get_uniref_mapping(gene_list):
    print("getting uniref mappings for {0}".format(gene_list))
    return ([])


def collect_uniref_sequences(uniref_genes, output_file):
    print("collecting uniref representative sequences in {0} into {1}".format(uniref_genes,
                                                                              output_file))


# This scans through all genes in the list model_genes and collects those protein
# sequences from the model_gene proteome

def collect_model_gene_sequences(model_genes, output_file):
    # print( "collecting model gene sequences in {0} into {1}".
    #          format( model_genes, output_file ) )

    # Its assumed here that it will be more efficient to retrieve all of the
    # relevant protein sequences from one genome in one call.  So we reorganize
    # the list accordingly:  First break apart the genome:gene names and then
    # collect in lists, by genome into a dict indexed by genome.

    genes_by_genome = {}
    for g in model_genes:
        genome, gene = g.split(':')
        if genes_by_genome.get(genome):
            genes_by_genome[genome] = genes_by_genome[genome] + [gene]
        else:
            genes_by_genome[genome] = [gene]

    # now do the calls, a genome at a time
    with open(output_file, "w") as pf:
        for geno in genes_by_genome:
            # print( "retrieving {0}".format( geno ) )
            sequences = fapi.get_protein_sequences(geno, genes_by_genome[geno])
            # print( sequences )
            pf.write(sequences)

    print("sequences written to {0}".format(output_file))


# presumbably protein sequences can be obtained from the query genome using
# the same api call

def get_query_genome_protein_sequences(genome_id, output_file):
    with open(output_file, "w") as pf:
        seqs = fapi.get_protein_sequences(genome_id)
        pf.write(seqs)

    print("query genome sequences written to {0}".format(output_file))


def get_min_index(a):
    min = a[0]
    min_i = 0
    for i in range(1, len(a)):
        if a[i] < min:
            min = a[i]
            min_i = i
    return (min_i)


# decided not to do this in an fapi since it doesn't touch the graphDB but
# maybe it should be

def find_best_homologs(query_seq_file, target_seq_file, noise_level=50, number_vals_to_report=5, threads=4):
    print("running blastp for {0} vs {1}".format(query_seq_file, target_seq_file))
    # wasn't able to get a pipe going
    # proc = subprocess.run( 'blastp -outfmt 6 -subject ' + target_seq_file +' -query ' + query_seq_file,
    #                       stdout=subprocess.PIPE, shell=True, universal_newline=True )
    #
    # for line in proc.stdout:
    #    print( line )

    # so using file output instead
    tmp_blast_output_file = os.path.join("blastp.results" + str(os.getpid()))

    blastp_cmd = f'blastp -outfmt 6 -subject {target_seq_file} -num_threads {threads} ' \
                 f'-query {query_seq_file} > {tmp_blast_output_file}'
    os.system(blastp_cmd)

    top_scores = [0] * number_vals_to_report
    top_records = [[]] * number_vals_to_report

    with open(tmp_blast_output_file) as bl:
        for line in bl:
            cols = line.strip().split()
            bl_score = float(cols[11])
            if bl_score > noise_level:
                min_i = get_min_index(top_scores)
                if bl_score > top_scores[min_i]:
                    top_scores[min_i] = bl_score
                    top_records[min_i] = cols

    return ([sorted(top_scores, reverse=True),
             sorted(top_records, reverse=True, key=lambda r: r[11])])


#
#  Here's the procedure to form the basis of the app method call.  Inputs are:
#
#    genome_id - some KBase-compliant id which can be used to retrieve protein sequences
#                 for now I'm using RefSeq ids ala 'GCF_000009625.1'
#
#    reaction_id - one of our reference reactions, i.e. 'rxn00371'
#
#    other_params - dict. of other inputs that could conceivably be additional app inputs
#                   such as similarity score thresholds and perhaps other filters TDB
#
def find_candidate_genes_for_reaction(reaction_id, genome_id, max_candidates=5):
    print("this is find_candidate_genes_for_reaction {0} in genome {1}".
          format(reaction_id, genome_id))

    # This is included here to show where I would do the reaction similarity.
    # Right now its set to return those reactions with either score == 1.0
    # If you want to avoid this, just remove the loop and do this
    #  reaction_model_genes = get_model_genes_for_reaction( reaction_id )

    reaction_set = get_similar_reactions(reaction_id)
    print("reaction set is {0}".format(reaction_set))
    reaction_model_genes = []
    for rxn in reaction_set:
        reaction_model_genes = reaction_model_genes + get_model_genes_for_reaction(rxn)

    # remove any duplicates
    reaction_model_genes = list(set(reaction_model_genes))

    using_uniref_mapping = False

    target_gene_seq_file = "targ_genes.fasta"

    if using_uniref_mapping:
        uniref_genes = get_uniref_mapping(reaction_model_genes)
        collect_uniref_sequences(uniref_genes, target_gene_seq_file)
    else:
        collect_model_gene_sequences(reaction_model_genes, target_gene_seq_file)

    query_genome_prot_seq_file = "query_genome_genes.fasta"
    get_query_genome_protein_sequences(genome_id, query_genome_prot_seq_file)

    results = find_best_homologs(query_genome_prot_seq_file, target_gene_seq_file,
                                 number_vals_to_report=max_candidates)
    return (results)


# Main program:  make one call for demo

if __name__ == "__main__":

    print("program starts")

    top_scores, top_hits = find_candidate_genes_for_reaction('rxn00371', 'GCF_000009625.1')

    print("results:")
    for h in top_hits:
        print(h)

    print("program ends")

