#!/usr/bin/env bash

set -u

DIR=$HOME/classes/INFO-492/rust-strobemers/data/experiment2
OUTDIR=$(mktemp -d)

for NUM in $(seq 0 9); do
    cargo run --bin kmer_comparison -- -k 5 \
        "$DIR/references${NUM}.fasta" \
        "$DIR/query${NUM}.fasta" \
        "$OUTDIR/out${NUM}.csv"
done

OUTFILE="$DIR/kmer.csv"
head -n1 $OUTDIR/out0.csv    > $OUTFILE
tail -n +2 -q $OUTDIR/out* >> $OUTFILE

rm -rf "$OUTDIR"

echo "Done, see outfile \"$OUTFILE\""