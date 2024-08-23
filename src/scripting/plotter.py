import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import os

import tools

def plot(experiment_name):
    relative_path_to_script = tools.experiment_dir(experiment_name)

    fontdict = {'fontsize': 15}
    plt.figure(figsize=(6, 6))

    # TODO
    estimation_files = glob.glob(os.path.join(relative_path_to_script, "query*.fasta"))
    alignment_files = glob.glob(os.path.join(relative_path_to_script, "references*.fasta"))
    for q, r in zip(estimation_files, alignment_files):
        data = pd.read_csv(
            os.path.join(relative_path_to_script, "estimation-output.csv")
        )
        alignment_data = pd.read_csv(
            os.path.join(relative_path_to_script, "alignment-output.csv")
        )
        data["edit_distance"] = alignment_data["edit_distance"]
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
