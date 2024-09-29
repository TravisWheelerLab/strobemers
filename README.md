# Fast alignment-free similarity estimation methods in Rust.

## Overview
This repository contains command-line executables that estimate the edit distance between query and reference sequences using:
1. Jaccard similarity of k-mer bags.
2. Jaccard similarity of strobemer bags.
3. maybe in the future, tensor slide sketch

## Progress
- Unit tests for everything but strobemers

## Development
- Develop unit tests for strobemers
- Write wrapper script to automate experiments
