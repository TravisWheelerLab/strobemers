import os

def run_strobemer(
    references_file,
    query_file,
    protocol,
    sim_method,
    order,
    strobe_length,
    strobe_window_gap,
    strobe_window_length,
    step
):
    command = (
        "cargo run --bin strobemer_comparison -- "
        f"{references_file} "
        f"{query_file} "
        f"-p {protocol} "
        f"-m {sim_method} "
        f"-o {order} "
        f"-l {strobe_length} "
        f"--w-gap {strobe_window_gap} "
        f"--w-len {strobe_window_length} "
        f"-s {step} "
    )
    print(f"Executing: {command}")
    os.system(command)

def strobemer():
    for order in [2]:
        for strobe_length in [5]:
            for strobe_window_gap in [0]:
                for strobe_window_length in [40]:
                    for step in [2, 3]:
                        run_strobemer(
    "artificial.fasta", "DF000000975.fa", "rand", "jaccard_similarity",
    order, strobe_length, strobe_window_gap, strobe_window_length, step)
strobemer()