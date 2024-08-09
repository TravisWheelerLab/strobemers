use std::time::Instant;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use bio::io::fasta;
use anyhow::Result;
use clap::Parser;

use alignment_free_methods;

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


fn create_csv_with_headers(name: String) -> Result<File>{
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let mut file = File::create(format!("{}/data/outputs/{}.csv", project_dir, name))?;
    writeln!(file, "ref_name,query_name,seed_name,estimation,ref_time,query_time")?;
    Ok(file)
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: KmerArgs) -> Result<()> {
    let mut csv_file = create_csv_with_headers(String::from("kmer-output"))?;
    let seed_name = format!("{}-mers", &args.k);
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = fasta::Reader::from_file(
        Path::new(&project_dir)
        .join("data/inputs")
        .join(&args.common.query_file)
    )?;

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

        let reference_reader = fasta::Reader::from_file(
            Path::new(&project_dir)
            .join("data/inputs")
            .join(&args.common.references_file)
        )?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            let reference_seed_start_time = Instant::now();
            let reference_seeds = alignment_free_methods::seq_to_kmers(
                reference_record.seq(),
                args.k.clone(),
                1
            )?;
            let reference_time = reference_seed_start_time.elapsed().as_secs_f64();

            let estimation: f64 = alignment_free_methods::jaccard_similarity(
                &reference_seeds,
                &query_seeds,
            )?;

            print!("\r comparisons done: {:?}", i);
            io::stdout().flush().unwrap();

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