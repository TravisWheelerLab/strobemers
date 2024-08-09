# Fast alignment-free similarity estimation methods in Rust.

## Overview
This repository contains command-line executables that estimate the edit distance between query and reference sequences using a few kinds of seeds. The executables take a reference fasta and a query fasta as positional arguments, and seed-specific flagged arguments. The following seed kinds have executables:
1. k-mers
2. spaced k-mers
3. strobemers

Similarity is calculated as the Jaccard similarity between seed-bags. Results are outputted to a directory specified by the "-e \[experiment-name\]" flag (specifically, data/outputs/\[experiment-name\]/\[seed-kind\]-output.csv).

I have not implemented subsampling techniques, and probably won't get to that by August 26.

## Progress and Development
* Check that bit-operations work.
* Generate data and create graphs.
