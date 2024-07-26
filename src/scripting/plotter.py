import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

def plot_line(data_file, method_label):
    data = pd.read_csv(data_file)
    plt.errorbar(x=data["edit distance"] / 1000, y=data["mean"],
        yerr=[data["mean"] - data["lower confidence bound"], data["upper confidence bound"] - data["mean"]],
        fmt='o',
        elinewidth=2,
        capsize=4,
        label=f'{method_label}')
    plt.grid(True)


def plot_log_line(data_file, method_label):
    data = pd.read_csv(data_file)
    plt.errorbar(x=data["edit distance"] / 1000, y=np.log(data["mean"]),
        yerr=[np.log(data["mean"]) - np.log(data["lower confidence bound"]), np.log(data["upper confidence bound"]) - np.log(data["mean"])],
        fmt='o',
        elinewidth=2,
        capsize=4,
        label=f'{method_label}')
    plt.grid(True)

def gapmer_line(db_size, k, gaps, step):
    data_file = ("../../tests/outputs/"
        f"data_{db_size}/"
        "gapmer/"
        "jaccard_similarity/"
        f"k_{k}/"
        f"gaps_{gaps}/"
        f"step_{step}/"
        "data.csv"
    )
    plot_line(data_file, f"{gaps}-gap {k}-mers, step={step}")

def alignment_line():
    data_file = "../../tests/outputs/" +\
        "data_2001/" +\
        "alignment/" +\
        "data.csv"

    plot_line(data_file, "Levenshtein edit distance")


def kmer_line(db_size, k, step):
    data_file = "../../tests/outputs/" +\
        f"data_{db_size}/" +\
        "kmer/" +\
        "jaccard_similarity/" +\
        f"k_{k}/" +\
        f"step_{step}/" +\
        "data.csv"
    plot_line(data_file, f"{k}-mers, step={step}")

def spaced_kmer_line(db_size, k, spaces, step):
    data_file = "../../tests/outputs/" +\
        f"data_{db_size}/" +\
        "spaced_kmer/" +\
        "jaccard_similarity/" +\
        f"k_{k}/" +\
        f"spaces_{spaces}/" +\
        f"step_{step}/" +\
        "data.csv"
    plot_line(data_file, f"{spaces}-spaced {k}-mers, step={step}")


def plot_kmers():
    for k in [5, 6, 7]:
        kmer_line(db_size=2001, k=k, step=1)

def plot_spaced_kmers():
    for k in [5, 6, 8]:
        for spaces in [1, 2]:
            spaced_kmer_line(2001, k, spaces, 1)
    
def strobemer_line(db_size, order, strobe_length, strobe_window_gap, strobe_window_length, step):
    data_file = ("../../tests/outputs/"
        f"data_{db_size}/"
        "strobemer/"
        "jaccard_similarity/"
        f"order_{order}/"
        f"strobe_length_{strobe_length}/"
        f"strobe_window_gap_{strobe_window_gap}/"
        f"strobe_window_length_{strobe_window_length}/"
        f"step_{step}/"
        "data.csv"
    )
    plot_log_line(data_file, f"({order}, {strobe_length}, {strobe_window_gap}, {strobe_window_length}, {step})-strobemer")


def plot_strobemers():
    for order in [2]:
        for strobe_length in [5]:
            for strobe_window_gap in [0]:
                for strobe_window_length in [40]:
                    for step in [1]:
                        strobemer_line(2001, order, strobe_length, strobe_window_gap, strobe_window_length, step)

def main():
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(4, 4))
    plt.title("Methods v. edit distance", fontdict=fontdict)
    plt.xlabel("True edit dissimilarity", fontdict=fontdict)
    plt.ylabel("Log(Estimation method)", fontdict=fontdict)

    # Put lines here
    alignment_line()
    plt.legend()
    plt.show()

    plt.legend()
    plt.show()

main()