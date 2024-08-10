import os
import os
def run_alignment(experiment_name, experiment_version):
    command = (
        "cargo run --bin alignment -- "
        f"references.fasta "
        f"query.fasta "
        f"-e {experiment_name} "
        f"-v {experiment_version}"
    )
    print(f"Executing: {command}")
    os.system(command)

run_alignment("experiment2", "_1")
