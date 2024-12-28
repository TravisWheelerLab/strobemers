#!/usr/bin/env bash

set -u

DIR=$HOME/classes/INFO-492/rust-strobemers/data/experiment2-ten-query
DATADIR=$DIR/data
SEQDIR=$DIR/sequences
FIGDIR=$DIR/figures
TEMPDIR=$(mktemp -d)

K=7
for NUM in $(seq 0 9); do
    cargo run --bin kmer_comparison -- -k $K \
        "$SEQDIR/references${NUM}.fasta" \
        "$SEQDIR/query${NUM}.fasta" \
        "$TEMPDIR/out${NUM}.csv"
done

OUTFILE="$DATADIR/${K}-mer.csv"
head -n1 $TEMPDIR/out0.csv    > $OUTFILE
tail -n +2 -q $TEMPDIR/out* >> $OUTFILE

rm -rf "$TEMPDIR"

echo "Done, see outfile \"$OUTFILE\""