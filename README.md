# Fast alignment-free similarity estimation methods in Rust.

## Overview
This repository contains command-line executables for estimating string similarity. Similarity is calculated as the Jaccard similarity of word-bag representations. Word types include:
1. k-mers
2. spaced k-mers
3. strobemers

Additionally, I may include subsampling techniques such as minimizers to speed up runtimes (I probably won't have time tho).

## Progress and Development
* Algorithms are mostly implemented (./src/lib.rs).
* Working on replicating the results of the SIM-R experiment in Kristoffer Sahlin's [paper](https://genome.cshlp.org/content/31/11/2080) and [code](https://github.com/ksahlin/strobemers). This will give me replicated data to work with.
