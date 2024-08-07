use std::time::Instant;
use std::io::{self, Write};
use std::path::Path;
use bio::io::fasta;
use bio::alignment::distance::levenshtein;
use anyhow::Result;
use clap::Parser;

use alignment_free_methods;

#[derive(Debug, Parser)]
struct StrobemerArgs {
    #[command(flatten)]
    common: alignment_free_methods::cli::CommonArgs,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(StrobemerArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: StrobemerArgs) -> Result<()> {
    let _seed_name = format!("alignment");
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = fasta::Reader::from_file(
        Path::new(&project_dir)
            .join("tests/inputs")
            .join(&args.common.query_file))?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let reference_reader = fasta::Reader::from_file(
            Path::new(&project_dir)
                .join("tests/inputs")
                .join(&args.common.references_file)
        )?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            print!("\rseeds generated:{:?}", i);
            io::stdout().flush().unwrap();
            
            let start = Instant::now();
            let _estimated_distance = levenshtein(
                query_record.seq(),
                reference_record.seq(),
            );
            let _duration = start.elapsed().subsec_millis();
        }
    }
    Ok(())
}