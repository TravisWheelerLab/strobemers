import random
import Levenshtein

""" This function creates two fastq files: a query fastq file with 1 random ATCG string length 1000,
and a reference fastq file with one modified query sequence every 10 edits from 0.0 dissimilarity
to 0.6 dissimilarity, or edit_distance in range(0, 600, 10).
The reference strings are generated from the query string with "random" edits.

This design creates a progression from similar to dissimilar strings based on actual edit distance,
not just number of edits. The only bias is the data is that there is only one query sequence, so I
recommend using multiple query and reference sequences in your experiments.
"""
def create_some_data(query_file, reference_file):
    directory = "data/inputs/simulated_data/"
    seq1 = "".join([random.choice("ACGT") for i in range(1000)])
    with open(directory + query_file, 'w') as o:
        o.write(f">query1 | len_1000 \n{seq1}")

    with open(directory + reference_file, 'w') as o:
        for edit_distance in range(0, 600, 10):
            seq2 = list(seq1)
            required_distance = edit_distance
            while required_distance > 0:
                idx = random.choice(range(0, len(seq2) - 1))
                mutation_type = random.choice(['deletion', 'substitution', 'insertion'])
                if mutation_type == 'deletion':
                    del seq2[idx]
                elif mutation_type == 'substitution':
                    seq2[idx] = random.choice("ACGT")
                elif mutation_type == 'insertion':
                    seq2.insert(idx, random.choice("ACGT"))
                current_distance = Levenshtein.distance(seq1, seq2)
                required_distance = edit_distance - current_distance
            o.write(f">query1+{edit_distance} | len_{len(seq2)} +{edit_distance} \n{''.join(seq2)}\n")
            print(f"string {edit_distance/10} generated", end='\r', flush=True)
    

create_some_data("query1.fasta", "references1.fasta")

