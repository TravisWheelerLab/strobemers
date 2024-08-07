use std::time::Instant;
use std::io::{self, Write};
use std::path::Path;
use bio::io::fastq;
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

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: KmerArgs) -> Result<()> {
    let _seed_name = format!("{}-mers", &args.k);
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = fastq::Reader::from_file(
        Path::new(&project_dir).join(&args.common.query_file)
    )?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seeds = alignment_free_methods::seq_to_kmers(
            query_record.seq(),
            args.k.clone(),
            1
        )?;
        let reference_reader = fastq::Reader::from_file(
            Path::new(&project_dir).join(&args.common.references_file)
        )?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;
            //let reference_seeds: Vec<Vec<char>> = match seeds_db_conn.prepare(
            //    "SELECT seed FROM seeds WHERE seq_name = ?1 AND seed_name = ?2")?
            //    .exists(params![reference_record.id(), seed_name])? {
            let reference_seeds = alignment_free_methods::seq_to_kmers(
                reference_record.seq(),
                args.k.clone(),
                1
            )?;

            print!("\rseeds generated:{:?}", i);
            io::stdout().flush().unwrap();
            
            let _start = Instant::now();
            let _estimated_distance: f64 = alignment_free_methods::jaccard_similarity(
                &reference_seeds,
                &query_seeds,
            )?;
        }
    }
    Ok(())
}