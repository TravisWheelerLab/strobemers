use std::time::Instant;
use std::io::{self, Write};
use std::path::Path;
use bio::io::fasta;
use anyhow::Result;
use clap::Parser;

use alignment_free_methods;

#[derive(Debug, Parser)]
struct StrobemerArgs {
    #[command(flatten)]
    common: alignment_free_methods::cli::CommonArgs,

    #[arg(short='p', value_name = "STRING", help="Strobemer generation protocol", default_value = "rand")]
    protocol: String,
    #[arg(short='o', value_name = "INT", help="Strobemer order/how many strobes", default_value_t = 2)]
    order: usize,
    #[arg(short='l', value_name = "INT", help="Strobe length")]
    strobe_length: usize,
    #[arg(long="w_min", value_name = "INT", help="w_min: window selection parameter")]
    w_min: usize,
    #[arg(long="w_max", value_name = "INT", help="w_max: window selection parameter")]
    w_max: usize,
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
    let _seed_name = format!("({},{},{},{})-{}strobemers",
        &args.order,
        &args.strobe_length,
        &args.w_min,
        &args.w_max,
        &args.protocol
    );
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;

    let query_reader = fasta::Reader::from_file(
        Path::new(&project_dir)
            .join("tests/inputs")
            .join(&args.common.query_file))?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seeds = alignment_free_methods::seq_to_randstrobemers(
            query_record.seq(),
            args.order.clone(),
            args.strobe_length.clone(),
            args.w_min.clone(),
            args.w_max.clone(),
            1
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
            let reference_seeds = alignment_free_methods::seq_to_randstrobemers(
                reference_record.seq(),
                args.order.clone(),
                args.strobe_length.clone(),
                args.w_min.clone(),
                args.w_max.clone(),
                1
            )?;

            print!("\rseeds generated:{:?}", i);
            io::stdout().flush().unwrap();
            
            let start = Instant::now();
            let _estimated_distance: f64 = alignment_free_methods::jaccard_similarity(
                &reference_seeds,
                &query_seeds,
            )?;
            let _duration = start.elapsed().subsec_millis();
        }
    }
    Ok(())
}