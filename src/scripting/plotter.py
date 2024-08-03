import statistics
import random
import sqlite3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

def plot_line(conn, seed_name):
    query = f"""
        SELECT 
            q1.query_name, 
            q1.reference_name, 
            q1.score AS alignment_score, 
            q2.score AS method_score
        FROM 
            comparisons q1
        JOIN 
            comparisons q2 
        ON 
            q1.query_name = q2.query_name AND q1.reference_name = q2.reference_name
        WHERE 
            q1.seed_name = 'alignment' AND q2.seed_name = '{seed_name}';"""
    data = pd.read_sql_query(query, conn)
    
    sns.scatterplot(x='alignment_score', y=f'method_score', data=data)
    plt.legend(['(Alignment Score, Strobemers Score)'], loc='upper left')
    plt.grid(True)
    print(seed_name)

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

def alignment_line(conn):
    # Execute the query and load the result into a pandas DataFrame    
    plot_line(conn, "alignment")


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
    
def strobemer_line(conn, order, strobe_length, strobe_window_gap, strobe_window_length, step, method):

    plot_line(conn, f"({order},{strobe_length},{strobe_window_gap},{strobe_window_length},{step})-{method}strobemers")


def plot_strobemers(conn):
    for order in [2]:
        for strobe_length in [5]:
            for strobe_window_gap in [0]:
                for strobe_window_length in [40]:
                    for step in [3]:
                        for method in ['rand']:
                            strobemer_line(conn, order, strobe_length, strobe_window_gap, strobe_window_length, step, method)

def main():
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(4, 4))
    plt.title("Methods v. edit distance", fontdict=fontdict)
    plt.xlabel("True edit dissimilarity", fontdict=fontdict)
    plt.ylabel("Log(Estimation method)", fontdict=fontdict)

    conn = sqlite3.connect('tests/outputs/comparisons.db')
    plot_strobemers(conn)
    plt.yscale('log')


    plt.legend()
    plt.show()

    plt.legend()
    plt.show()

main()