def run_spaced_kmer(
    db_size,
    sim_method,
    k,
    spaces,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"spaced_kmer/{sim_method}/"
        f"k_{k}/"
        f"spaces_{spaces}/"
        f"step_{step}/"
    )
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method spaced_kmer "
        f"--similarity-method {sim_method} "
        f"-k {k} "
        f"--spaces {spaces} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def spaced_kmer():
    for k in [5, 6]:
        for spaces in [2]:
            run_spaced_kmer(2001, "jaccard_similarity", k, spaces, 1)

