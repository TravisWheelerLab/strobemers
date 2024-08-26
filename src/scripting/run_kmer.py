import os

def kmers():
    run_kmer(
        references_file="simulated_data/references1.fasta",
        query_file="simulated_data/query1.fasta",
        k=30,
    )

kmers()
