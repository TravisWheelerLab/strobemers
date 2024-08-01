use std::time::Instant;
use std::io::{self, Write};
use std::path::Path;
use bio::io::fasta;
use anyhow::Result;
use clap::Parser;
use rusqlite::{self, params};

use alignment_free_methods;

#[derive(Debug, Parser)]
struct StrobemerArgs {
    #[command(flatten)]
    common: alignment_free_methods::cli::CommonArgs,

    #[arg(short='p', value_name = "STRING")]
    protocol: String,
    #[arg(short='o', value_name = "INT")]
    order: usize,
    #[arg(short='l', value_name = "INT", long)]
    strobe_length: usize,
    #[arg(long="w-gap", value_name = "INT")]
    strobe_window_gap: usize,
    #[arg(long="w-len", value_name = "INT")]
    strobe_window_length: usize,
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
    let seed_name = format!("({},{},{},{},{})-{}strobemers",
        &args.order,
        &args.strobe_length,
        &args.strobe_window_gap,
        &args.strobe_window_length,
        &args.common.step,
        &args.protocol
    );
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let comparison_db_conn = alignment_free_methods::cli::initialize_comparison_db(
        Path::new(&project_dir).join("tests/outputs/comparisons.db")
    )?;

    let query_reader = fasta::Reader::from_file(
        Path::new(&project_dir)
            .join("tests/inputs")
            .join(&args.common.query_file))?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seeds = alignment_free_methods::generate_strobemers(
            query_record.seq(),
            args.order.clone(),
            args.strobe_length.clone(),
            args.strobe_window_gap.clone(),
            args.strobe_window_length.clone(),
            args.common.step.clone(),
            None
        )?;
        let reference_reader = fasta::Reader::from_file(
            Path::new(&project_dir)
                .join("tests/inputs")
                .join(&args.common.references_file)
        )?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;
            //let reference_seeds: Vec<Vec<char>> = match seeds_db_conn.prepare(
            //    "SELECT seed FROM seeds WHERE seq_name = ?1 AND seed_name = ?2")?
            //    .exists(params![reference_record.id(), seed_name])? {
            let reference_seeds: Vec<Vec<char>> = alignment_free_methods::generate_strobemers(
                reference_record.seq(),
                args.order.clone(),
                args.strobe_length.clone(),
                args.strobe_window_gap.clone(),
                args.strobe_window_length.clone(),
                args.common.step.clone(),
                None
            )?;

            print!("\rseeds generated:{:?}", i);
            io::stdout().flush().unwrap();
            
            let start = Instant::now();
            let estimated_distance: f64 = alignment_free_methods::jaccard_similarity(
                &reference_seeds,
                &query_seeds,
            )?;
            let duration = start.elapsed().subsec_millis();
            comparison_db_conn.execute(
                "INSERT OR REPLACE INTO comparisons (query_name, reference_name, seed_name, score, time) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![query_record.id(), reference_record.id(), seed_name, estimated_distance, duration],
            )?;
        }
    }
    Ok(())
}