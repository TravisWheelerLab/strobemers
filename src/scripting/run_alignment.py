import os
import os
def run_alignment(references_file, query_file):

    command = (
        "cargo run --bin alignment -- "
        f"{references_file} "
        f"{query_file} "
    )
    print(f"Executing: {command}")
    os.system(command)

def alignment():
    run_alignment("simulated_data/references1.fasta", "simulated_data/query1.fasta")
    
alignment()
