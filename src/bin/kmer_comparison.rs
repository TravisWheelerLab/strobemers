use std::time::Instant;
use anyhow::Result;
use clap::Parser;

use alignment_free_methods::cli::*;
use alignment_free_methods::{seq_to_kmers, jaccard_similarity};

#[derive(Debug, Parser)]
struct KmerArgs {
    #[command(flatten)]
    common: alignment_free_methods::cli::CommonArgs,

    #[arg(short='k', value_name = "INT", help="Substring length")]
    k: usize,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(KmerArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run(args: KmerArgs) -> Result<()> {
    let seed_name = format!("{}-mers", &args.k);
    let query_reader = create_reader(&args.common.query_file)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        // There should only be 1!
        let query_record = query_record?;
        let query_seed_start_time = Instant::now();
        let query_seeds = alignment_free_methods::seq_to_kmers(
            query_record.seq(),
            args.k.clone(),
            1
        )?;
        let query_time = query_seed_start_time.elapsed().as_secs_f64();

        let reference_reader = create_reader(&args.common.references_file)?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            let reference_seed_start_time = Instant::now();
            let reference_seeds = seq_to_kmers(
                reference_record.seq(),
                args.k.clone(),
                1
            )?;

            let reference_time = reference_seed_start_time.elapsed().as_secs_f64();
            let estimation = jaccard_similarity(&reference_seeds,&query_seeds,)?;
            print_update(i);
        }
    }
    Ok(())
}