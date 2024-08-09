use std::time::Instant;
use std::fs::File;
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

fn create_csv_with_headers(name: String) -> Result<File>{
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let mut file = File::create(format!("{}/data/outputs/{}.csv", project_dir, name))?;
    writeln!(file, "ref_name,query_name,seed_name,edit_distance,edit_distance_time")?;
    Ok(file)
}
// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: StrobemerArgs) -> Result<()> {
    let mut csv_file = create_csv_with_headers(String::from("alignment-output"))?;
    let seed_name = format!("alignment");
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = fasta::Reader::from_file(
        Path::new(&project_dir)
        .join("data/inputs")
        .join(&args.common.query_file)
    )?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let reference_reader = fasta::Reader::from_file(
            Path::new(&project_dir)
                .join("data/inputs")
                .join(&args.common.references_file)
        )?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            print!("\rseeds generated:{:?}", i);
            io::stdout().flush().unwrap();

            let start = Instant::now();
            let estimated_distance = levenshtein(
                query_record.seq(),
                reference_record.seq(),
            );
            let duration = start.elapsed().as_secs_f64();

            writeln!(csv_file, "{},{},{},{},{}",
                reference_record.id(),
                query_record.id(),
                seed_name,
                estimated_distance,
                duration
            )?;
        }
    }
    Ok(())
}