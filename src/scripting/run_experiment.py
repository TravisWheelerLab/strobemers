import pathlib
import os

import tools
import plotter
import run_strobemers
import run_alignment
import simulate_data

def fasta_file_paths(experiment_name):
    relative_path_to_script = tools.experiment_dir(experiment_name)
    query = os.path.join(relative_path_to_script, "query.fasta")
    references = os.path.join(relative_path_to_script, "references.fasta")
    return (query, references)


# Create a relative path from the run directory to the script
def experiment1():
    experiment_name = "experiment1"
    query, references = fasta_file_paths(experiment_name)
    os.mkdir(os.path.dirname(query))

    simulate_data.create_some_data(experiment_name, 1000, range(0, 500, 10))
    run_strobemers.run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=7,
        w_min=1, w_max=25,
        experiment_name=experiment_name
    )
    run_alignment.run_alignment(experiment_name)

    plotter.plot(experiment_name)

# Create a relative path from the run directory to the script
def experiment2():
    experiment_name = "experiment2"
    query, references = fasta_file_paths(experiment_name)
    os.mkdir(os.path.dirname(query))

    simulate_data.create_correlated_data(experiment_name, 2, 10000, range(0, 5000, 10))

    run_strobemers.run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=7,
        w_min=1, w_max=25,
        experiment_name=experiment_name
    )
    run_alignment.run_alignment(experiment_name)

    plotter.plot(experiment_name)

experiment2()