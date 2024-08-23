import os
import glob

import tools

def run_strobemer(
    protocol,
    order,
    strobe_length,
    w_min,
    w_max,
    experiment_name
):
    relative_path_to_script = tools.experiment_dir(experiment_name)
    query_files = glob.glob(os.path.join(relative_path_to_script, "query*.fasta"))
    references_files = glob.glob(os.path.join(relative_path_to_script, "references*.fasta"))
    for q, r in zip(sorted(query_files), sorted(references_files)):
        command = (
            "cargo run --bin strobemer_comparison -- "
            f"{os.path.basename(r)} "
            f"{os.path.basename(q)} "
            f"-p {protocol} "
            f"-o {order} "
            f"-l {strobe_length} "
            f"--w-min {w_min} "
            f"--w-max {w_max} "
            f"-e {experiment_name} "
        )
        print(f"Executing: {command}")
        os.system(command)
