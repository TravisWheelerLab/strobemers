
def run_kmer(
    db_size,
    sim_method,
    k,
    step
):
    outpath = (f"../../tests/outputs/"
        f"data_{db_size}/"
        f"kmer/{sim_method}/"
        f"k_{k}/"
        f"step_{step}/"
    )
    print(f"Creating output directory...")
    os.system(f"mkdir -p {outpath}")
    command = (
        "cargo run --bin generate_comparison_data -- "
        f"-i ../../tests/inputs/sequences_{db_size}.csv "
        f"-o {outpath}data.csv "
        f"--representation-method kmer "
        f"--similarity-method {sim_method} "
        f"-k {k} "
        f"--step {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def kmers(k, step):
    run_kmer(
        db_size=2001,
        sim_method="jaccard_similarity",
        k=k,
        step=step
    )

