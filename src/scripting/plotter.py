import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

# takes a pandas dataframe with "estimation" and "edit_distance" columns
# and draws the points on a pre-existing plot.
def plot_line(data):
    plt.scatter(x='edit_distance', y='estimation', data=data)
    plt.grid(True)

def plot_kmers():
    estimation_data = pd.read_csv(f"data/outputs/{experiment_number}/kmer-output.csv")
    alignment_data = pd.read_csv(f"data/outputs/{experiment_number}/alignment-output.csv")
    estimation_data["edit_distance"] = alignment_data["edit_distance"]
    plot_line(estimation_data)

def plot_strobemers(experiment_name):
    estimation_data = pd.read_csv(f"data/outputs/{experiment_name}/strobemer-output.csv")
    alignment_data = pd.read_csv(f"data/outputs/{experiment_name}/alignment-output.csv")
    estimation_data["edit_distance"] = alignment_data["edit_distance"]
    plot_line(estimation_data)

def plot(experiment_name, seed_name):
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(6, 6))
    plt.title(f"{seed_name} similarity vs. edit distance", fontdict=fontdict)
    plt.xlabel("True edit distance", fontdict=fontdict)
    plt.ylabel(f"Log-scaled jaccard-similarity of \n {seed_name} bags", fontdict=fontdict)

    plot_strobemers(experiment_name)
    plt.yscale('log')

    plt.legend()
    plt.show()

def experiment1():
    plot("experiment1", "randstrobes")

def experiment1_1():
    plot("experiment1_1", "randstrobes")

def experiment1_2():
    plot("experiment1_2", "randstrobes")

def experiment1_3():
    plot("experiment1_3", "(2,15,0,30)-minstrobe")

def experiment2_1():
    plot("experiment2_1", "(2,15,0,30)-randstrobe")

def experiment2_2():
    plot("experiment2_2", "randstrobes")

def experiment2_3():
    plot("experiment2_3", "randstrobes")

def experiment2_4():
    plot("experiment2_4", "randstrobes")


def experiment3_1():
    plot("experiment3_1", "randstrobes")

experiment1_1()
