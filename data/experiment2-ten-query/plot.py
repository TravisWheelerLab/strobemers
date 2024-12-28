#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np

SUFFICIENTLY_SMALL = 5e-4

def preprocess(df):
    df['edit_distance'] = df['ref_name'].apply(lambda x: int(re.search(r'\+(\d+)', x).group(1)))
    df['log_estimation'] = np.log(df['estimation'].apply(lambda x: x if x else SUFFICIENTLY_SMALL))
    df['is_zero'] = df['log_estimation'] == np.log(SUFFICIENTLY_SMALL)


DIR = "/home/zoc/classes/INFO-492/rust-strobemers/data/experiment2-ten-query"


#csv_file = DIR + "/data/(2,5,5,25)-minstrobemer.csv"
csv_file = DIR + "/data/7-mer.csv"
df = pd.read_csv(csv_file)

preprocess(df)

plt.figure(figsize=(8, 6))
plt.scatter(df['edit_distance'][~df['is_zero']], df['log_estimation'][~df['is_zero']], color='blue', s=6)
plt.scatter(df['edit_distance'][df['is_zero']], df['log_estimation'][df['is_zero']], color='red', s=6)

plt.xlabel('Edit Distance')
plt.ylabel('log_e(Estimation)')

N_QUERIES = 10
SEQ_BASE_LEN = 1000 # 
plt.title(f"Edit Distance vs {df['seed_name'][0]} Estimation ({N_QUERIES}x{SEQ_BASE_LEN}-bp queries)")

fig_name = f"{DIR + '/figures/' + df['seed_name'][0]}.png"
plt.savefig(fig_name, format='png')