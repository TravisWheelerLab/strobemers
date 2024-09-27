import pathlib
import os
import shutil
import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import glob

import tools
import simulate_data

def fasta_file_paths(experiment_name):
    relative_path_to_script = tools.experiment_dir(experiment_name)
    query = os.path.join(relative_path_to_script, "query.fasta")
    references = os.path.join(relative_path_to_script, "references.fasta")
    return (query, references)

def fasta_file_lists(experiment_name):
    relative_path_to_script = tools.experiment_dir(experiment_name)
    query_files = glob.glob(os.path.join(relative_path_to_script, "query*.fasta"))
    references_files = glob.glob(os.path.join(relative_path_to_script, "references*.fasta"))
    return (query_files, references_files)

def run_alignment(experiment_name):
    query_files, references_files = fasta_file_lists(experiment_name)
    for q, r in zip(sorted(query_files), sorted(references_files)):
        command = (
            "cargo run --bin alignment -- "
            f"{os.path.basename(r)} "
            f"{os.path.basename(q)} "
            f"-e {experiment_name} "
        )
        print(f"Executing: {command}")
        os.system(command)


def run_kmer(k, experiment_name):
    query_files, references_files = fasta_file_lists(experiment_name)
    for q, r in zip(sorted(query_files), sorted(references_files)):
        command = (
            "cargo run --bin kmer_comparison -- "
            f"{os.path.basename(r)} "
            f"{os.path.basename(q)} "
            f"-k {k} "
            f"-e {experiment_name} "
        )
        print(f"Executing: {command}")
        os.system(command)

def run_strobemer(
    protocol,
    order,
    strobe_length,
    w_min,
    w_max,
    experiment_name,
):
    query_files, references_files = fasta_file_lists(experiment_name)
    for i, (q, r) in enumerate(zip(sorted(query_files), sorted(references_files))):
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
            f"--run-alignment "
        )
        print(f"Executing: {command}")
        os.system(command)


def plot(experiment_name):
    relative_path_to_script = tools.experiment_dir(experiment_name)

    fontdict = {'fontsize': 15}
    plt.figure(figsize=(6, 6))

    data_files = glob.glob(os.path.join(relative_path_to_script, "output*.csv"))
    for data_file in data_files:
        data = pd.read_csv(
            os.path.join(relative_path_to_script, os.path.basename(data_file))
        )
        seed_name = data.at[0, "seed_name"]
        plt.scatter(x='edit_distance', y='estimation', data=data)

    plt.title(f"{seed_name} similarity vs. edit distance", fontdict=fontdict)
    plt.xlabel("True edit distance", fontdict=fontdict)
    plt.ylabel(f"Log-scaled jaccard-similarity of \n {seed_name} bags", fontdict=fontdict)
    plt.grid(True)
    if "strobe" in seed_name:
        plt.yscale('log')

    plt.legend()
    plt.savefig(os.path.join(relative_path_to_script, "figure.png"))

    plt.show()

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
    shutil.rmtree(os.path.dirname(query))
    os.mkdir(os.path.dirname(query))

    n_pairs = 10
    simulate_data.create_correlated_data(experiment_name, n_pairs, 1000, range(0, 500, 10))

    relative_path_to_script = tools.experiment_dir(experiment_name)
    string_pairs = glob.glob(os.path.join(relative_path_to_script, "query*.fasta"))
    run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=7,
        w_min=1, w_max=25,
        experiment_name=experiment_name
    )
    run_alignment(experiment_name)

    plot(experiment_name)

experiment2()