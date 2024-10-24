# Fast alignment-free similarity estimation methods in Rust.

## Overview
This repository contains command-line executables that estimate the edit distance between query and reference sequences using:
1. Jaccard similarity of k-mer bags.
2. Jaccard similarity of strobemer bags.

## Progress
- Minstrobes implemented and tested.
- 10,000-nt sequences with pre-computed edit distance generated for comparison with seed-bag-comparison estimates ((/data/experiment3)[`data/experiment3_`]).

## Development
- Randstrobes and Hybridstrobes
