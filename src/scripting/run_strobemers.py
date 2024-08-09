import os

def run_strobemer(
    references_file,
    query_file,
    protocol,
    order,
    strobe_length,
    w_min,
    w_max,
    experiment_name
):
    command = (
        "cargo run --bin strobemer_comparison -- "
        f"{references_file} "
        f"{query_file} "
        f"-p {protocol} "
        f"-o {order} "
        f"-l {strobe_length} "
        f"--w-min {w_min} "
        f"--w-max {w_max} "
        f"-e {experiment_name}"
    )
    print(f"Executing: {command}")
    os.system(command)

def strobemer():
    for order in [2]:
        for strobe_length in [15]:
            for w_min in [1]:
                for w_max in [50]:
                    run_strobemer(
                        references_file="simulated_data/query1.fasta",
                        query_file="simulated_data/references1.fasta",
                        protocol="rand",
                        order=order,
                        strobe_length=strobe_length,
                        w_min=w_min, w_max=w_max,
                        experiment_name="experiment1"
    )
strobemer()