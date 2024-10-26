#!/usr/bin/env bash

set -u

DIR=$HOME/classes/INFO-498H/rust-strobemers/data/experiment3_100_query_higher_res
SEQDIR=$DIR/sequences
DATADIR=$DIR/data
FIGDIR=$DIR/figures
TEMPDIR=$(mktemp -d)

P="min"
O=2
L=6
WMIN=5
WMAX=25

for NUM in $(seq 0 9); do
    cargo run -r --bin strobemer_comparison -- \
        -p $P \
        -o $O  \
        -l $L \
        --w-min $WMIN \
        --w-max $WMAX \
        "$SEQDIR/references${NUM}.fasta" \
        "$SEQDIR/query${NUM}.fasta" \
        "$TEMPDIR/out${NUM}.csv"
done

OUTFILE="$DATADIR/($O,$L,$WMIN,$WMAX)-${P}strobemer.csv"
head -n1 $TEMPDIR/out0.csv    > $OUTFILE
tail -n +2 -q $TEMPDIR/out* >> $OUTFILE

rm -rf "$TEMPDIR"

echo "Done, see outfile \"$OUTFILE\""