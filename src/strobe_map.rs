use anyhow::{bail, Result};
use clap::Parser;
use csv::{ReaderBuilder, WriterBuilder};
use serde::Deserialize;
use std::collections::HashMap;
use std::io::{self,Write};

#[derive(Debug, Parser)]
#[command(about, author, version)]
/// Similarity Methods
struct Args {
    #[arg(short, long, help = "Path to a CSV of <string, string, integer>'s")]
    input: String,
    #[arg(short, long, help = "Path to an output CSV.", default_value = "tests/outputs/unnamed_data.csv")]
    output: String,
    #[arg(short, long, help = "Which string representation method (e.g. k-mers")]
    method: String,
    #[arg(short, long, help = "Distance estimation function (e.g. Jaccard similarity)")]
    distance_function: String,
    #[arg(long, help = "Spacing between representation windows", default_value_t = 1)]
    step: usize,
    #[arg(short, help = "For  k-mer-based methods")]
    k: Option<usize>,
    #[arg(long, help = "For minimizer-based methods")]
    minimizer_window_length: Option<usize>,
    #[arg(long, help = "For strobemer-based methods")]
    strobemer_order: Option<usize>,
    #[arg(long, help = "For strobemer-based methods")]
    strobe_length: Option<usize>,
    #[arg(long, help = "For strobemer-based methods")]
    strobe_window_gap: Option<usize>,
    #[arg(long, help = "For strobemer-based methods")]
    strobe_window_length: Option<usize>,
}
use similarity_methods;

/* Using serde to parse CSV data. */
#[derive(Debug, Deserialize)]
struct DatabaseRecord {
    base_sequence: String,
    modified_sequence: String,
    edit_distance: usize,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: Args) -> Result<()> {
    // Rolling magnitude and variability data.
    let mut edit_distance_sums: HashMap<usize, f64> = HashMap::new();
    let mut edit_distance_squared_sums: HashMap<usize, f64> = HashMap::new();
    let mut edit_distance_counts: HashMap<usize, f64> = HashMap::new();
    let mut rdr = ReaderBuilder::new().from_path(&args.input)?;

    for (i, result) in rdr.deserialize().enumerate() {
        let record: DatabaseRecord = result?;
        let base_seq: Vec<char> = record.base_sequence.chars().collect();
        let mod_seq: Vec<char> = record.modified_sequence.chars().collect();        
        let estimated_distance: f64 = match args.method.as_str() {
            "kmer" => similarity_methods::kmer_similarity(
                &base_seq,
                &mod_seq,
                &args.distance_function,
                args.k.clone()
                    .expect("argument 'k' not provided!"),
                args.step.clone()
            )?,
            "minimizer" => similarity_methods::minimizer_similarity(
                &base_seq,
                &mod_seq,
                &args.distance_function,
                args.k.clone()
                    .expect("argument 'k' not provided!"),
                args.minimizer_window_length.clone()
                    .expect("argument 'minimizer_window_length' not provided!"),
                args.step.clone()
            )?,
            "strobemer" => similarity_methods::strobemer_similarity(
                &base_seq,
                &mod_seq,
                &args.distance_function,
                args.strobemer_order.clone()
                    .expect("argument 'strobemer_order' not provided!"),
                args.strobe_length.clone()
                    .expect("argument 'strobe_length' not provided!"),
                args.strobe_window_gap.clone()
                    .expect("argument 'strobe_window_gap' not provided!"),
                args.strobe_window_length.clone()
                    .expect("argument 'strobe_window_length' not provided!"),
                args.step.clone()
            )?,
            _ => {
                bail!("Unknown representation method: {}", args.method.as_str());
            }
        };

        let sum = edit_distance_sums.entry(record.edit_distance.clone()).or_insert(0.0);
        let squared_sum = edit_distance_squared_sums.entry(record.edit_distance.clone()).or_insert(0.0);
        let count = edit_distance_counts.entry(record.edit_distance.clone()).or_insert(0.0);
        *sum += estimated_distance;
        *squared_sum += estimated_distance * estimated_distance;
        *count += 1.0;

        print!("\rString pair #{:?} has been processed!", i);
        io::stdout().flush().unwrap();
    }

    // Compute mean and confidence intervals from rolling data.
    let mut lower_bounds: HashMap<usize, f64> = HashMap::new();
    let mut upper_bounds: HashMap<usize, f64> = HashMap::new();
    let mut means: HashMap<usize, f64> = HashMap::new();
    for (edit_distance, count) in &edit_distance_counts {
        let sum = edit_distance_sums.get(&edit_distance).unwrap();
        let squared_sum = edit_distance_squared_sums.get(&edit_distance).unwrap();

        let mean = sum / count;
        let variance = (squared_sum - mean * sum) / count;
        let mean_se = (variance / count).sqrt();

        means.insert(edit_distance.clone(), mean);
        lower_bounds.insert(edit_distance.clone(), mean - mean_se);
        upper_bounds.insert(edit_distance.clone(), mean + mean_se);
    }

    // Create the header row for comparison data file
    let mut header_row = vec!["edit distance".to_string()];
    header_row.push("mean".to_string());
    header_row.push("lower confidence bound".to_string());
    header_row.push("upper confidence bound".to_string());
    header_row.push("number of string pairs".to_string());

    let mut wtr = WriterBuilder::new()
        .delimiter(b',')
        .from_path(&args.output)?;
    wtr.write_record(&header_row)?;

    let mut sorted_edit_distances: Vec<&usize> = edit_distance_counts.keys().collect();
    sorted_edit_distances.sort();
    for dist in sorted_edit_distances {
        let mut row = vec![dist.to_string()];
        row.push(means.get(&dist).unwrap().to_string());
        row.push(lower_bounds.get(&dist).unwrap().to_string());
        row.push(upper_bounds.get(&dist).unwrap().to_string());
        row.push(edit_distance_counts.get(&dist).unwrap().to_string());
        wtr.write_record(&row)?;
    }

    println!(r#"Done, see output "{}""#, args.output);
    Ok(())
}