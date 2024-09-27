use std::time::Instant;
use bio::alignment::distance::levenshtein;
use anyhow::Result;
use clap::Parser;

use alignment_free_methods::cli::*;

// --------------------------------------------------
fn main() {
    if let Err(e) = run(CommonArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: CommonArgs) -> Result<()> {
    let seed_name = format!("alignment");
    let query_reader = create_reader(&args.query_file)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let reference_reader = create_reader(&args.references_file)?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            print_update(i);
            let start = Instant::now();
            let estimation = levenshtein(query_record.seq(),reference_record.seq());
            let duration = start.elapsed().as_secs_f64();            
        }
    }
    Ok(())
}