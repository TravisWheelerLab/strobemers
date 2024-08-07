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
        references_file="data/sahlin_data/cDNA/ERR3588905_1.fastq",
        query_file="data/sahlin_data/cDNA/ERR3588905_1_full.fastq",
        k=30,
    )

kmers()
