import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

def plot_line(data):
    sns.scatterplot(x='edit_distance', y='estimation', data=data)
    plt.legend(['(Alignment Score, Strobemers Score)'], loc='upper left')
    plt.grid(True)


def plot_kmers():
    estimation_data = pd.read_csv(f"data/outputs/{experiment_number}/kmer-output.csv")
    alignment_data = pd.read_csv(f"data/outputs/{experiment_number}/alignment-output.csv")
    estimation_data["edit_distance"] = alignment_data["edit_distance"]
    plot_line(estimation_data)

def plot_strobemers(experiment_number):
    estimation_data = pd.read_csv(f"data/outputs/{experiment_number}/strobemer-output.csv")
    alignment_data = pd.read_csv(f"data/outputs/{experiment_number}/alignment-output.csv")
    print(estimation_data.head())
    print(alignment_data.head())
    edit_distance = alignment_data["edit_distance"]
    estimation_data["edit_distance"] = edit_distance
    plot_line(estimation_data)

def main():
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(4, 4))
    plt.title("Methods v. edit distance", fontdict=fontdict)
    plt.xlabel("True edit dissimilarity", fontdict=fontdict)
    plt.ylabel("Log(Estimation method)", fontdict=fontdict)

    plot_strobemers("experiment1")
    plt.yscale('log')

    plt.legend()
    plt.show()

main()