// WARNING: this binary does not implement functionality for the run_alignment flag.

use std::time::Instant;
use std::io::Write;
use alignment_free_methods::cli::create_query_reader;
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
    let mut csv_file = create_csv_with_estimation_headers(&args.common)?;
    let seed_name = format!("{}-mers", &args.k);
    let query_reader = create_query_reader(&args.common)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seed_start_time = Instant::now();
        let query_seeds = alignment_free_methods::seq_to_kmers(
            query_record.seq(),
            args.k.clone(),
            1
        )?;
        let query_time = query_seed_start_time.elapsed().as_secs_f64();

        let reference_reader = create_reference_reader(&args.common)?;
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
            writeln!(csv_file, "{},{},{},{},{},{:?}",
                reference_record.id(),
                query_record.id(),
                seed_name,
                estimation,
                reference_time,
                query_time
            )?;
        }
    }
    Ok(())
}