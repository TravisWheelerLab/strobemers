import random
import Levenshtein
import os

import tools

""" This function creates two fastq files: a query fastq file with 1 random ATCG string length 1000,
and a reference fastq file with one modified query sequence every 10 edits from 0.0 dissimilarity
to 0.6 dissimilarity, or edit_distance in range(0, 600, 10).
The reference strings are generated from the query string with "random" edits.

This design creates a progression from similar to dissimilar strings based on actual edit distance,
not just number of edits. The only bias is the data is that there is only one query sequence, so I
recommend using multiple query and reference sequences in your experiments.
"""
def create_some_data(experiment_name, seq_len, edit_distance_range):
    relative_path_to_script = tools.experiment_dir(experiment_name)
    query_file = os.path.join(relative_path_to_script, "query.fasta")
    reference_file = os.path.join(relative_path_to_script, "references.fasta")

    seq1 = "".join([random.choice("ACGT") for i in range(seq_len)])

    with open(query_file, 'w') as o:
        o.write(f">query1 | len_{seq_len} \n{seq1}")

    with open(reference_file, 'w') as o:
        for edit_distance in edit_distance_range:
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
    


def create_correlated_data(experiment_name, n, seq_len, edit_distance_range):
    relative_path_to_script = tools.experiment_dir(experiment_name)
    for n in range(n):
        query_file = os.path.join(relative_path_to_script, f"query{n}.fasta")
        references_file = os.path.join(relative_path_to_script, f"references{n}.fasta")


        seq1 = "".join([random.choice("ACGT") for i in range(seq_len)])
        with open(query_file, 'w') as o:
            o.write(f">query{n} | len_{seq_len} \n{seq1}")

        with open(references_file, 'w') as o:
            previous_distance = 0
            current_distance = 0
            seq2 = list(seq1)
            for edit_distance in edit_distance_range:
                required_distance = edit_distance - current_distance
                while required_distance > 0:
                    for _ in range(required_distance):
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

                o.write(f">query{n}+{current_distance} | len_{len(seq2)} +{edit_distance} \n{''.join(seq2)}\n")
                print(f"string {edit_distance} generated", end='\r', flush=True)