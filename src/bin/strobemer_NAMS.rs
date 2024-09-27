// WARNING: this binary does not implement functionality for the run_alignment flag.

use std::time::Instant;
use anyhow::Result;
use clap::Parser;

use alignment_free_methods::cli::*;
use alignment_free_methods::{seq_to_strobemers, jaccard_similarity};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(StrobemerArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run(args: StrobemerArgs) -> Result<()> {
    let seed_name = format_strobemer_seed_name(&args);

    let query_reader = create_reader(&args.common.query_file)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seed_start_time = Instant::now();
        let query_seeds = seq_to_strobemers(
            query_record.seq(),
            &args
        )?;
        let query_time = query_seed_start_time.elapsed().as_secs_f64();

        let reference_reader = create_reader(&args.common.references_file)?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            let reference_seed_start_time = Instant::now();
            let reference_seeds = seq_to_strobemers(
                reference_record.seq(),
                &args
            )?;
            let reference_time = reference_seed_start_time.elapsed().as_secs_f64();
            print_update(i);
            let estimation: f64 = jaccard_similarity(&reference_seeds, &query_seeds)?;
        }
    }
    Ok(())
}