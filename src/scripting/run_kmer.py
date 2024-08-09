import os
def run_kmer(references_file, query_file, k):

    command = (
        "cargo run --bin kmer_comparison -- "
        f"{references_file} "
        f"{query_file} "
        f"-k {k} "
    )
    print(f"Executing: {command}")
    os.system(command)

def kmers():
    run_kmer(
        references_file="simulated_data/references1.fasta",
        query_file="simulated_data/query1.fasta",
        k=30,
    )

kmers()
