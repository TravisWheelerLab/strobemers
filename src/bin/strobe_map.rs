use anyhow::Result;
use clap::{builder::PossibleValue, Parser, ValueEnum};
use csv::{ReaderBuilder, WriterBuilder};
use serde::Deserialize;
use similarity_methods::{
    kmer_similarity, minimizer_similarity, strobemer_similarity,
    DistanceFunction,
};
use std::{
    collections::HashMap,
    io::{self, Write},
};

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[arg(long, value_name = "METHOD")]
    method: SimilarityMethod,

    #[arg(short, long, value_name = "INPUT", help = "Input CSV file")]
    input: String,

    #[arg(
        short,
        long,
        default_value = "out.csv",
        value_name = "OUTPUT",
        help = "Output CSV file"
    )]
    output: String,

    #[arg(short, long, help = "Distance estimation function")]
    distance_function: DistanceFunction,

    #[arg(
        short,
        long,
        help = "Spacing between representation windows",
        default_value_t = 1
    )]
    step: usize,

    /// Size of k
    #[arg(short)]
    k: usize,

    /// Minimizer window length
    #[arg(long)]
    minimizer_window_length: usize,

    /// Strobemer Order
    #[arg(long)]
    strobemer_order: usize,

    /// Strobemer length
    #[arg(long)]
    strobemer_length: usize,

    /// Strobemer window gap
    #[arg(long)]
    strobemer_window_gap: usize,

    /// Strobemer window length
    #[arg(long)]
    strobemer_window_length: usize,
}

#[derive(Debug, Eq, PartialEq, Clone)]
enum SimilarityMethod {
    Kmer,
    Minimizer,
    Strobemer,
}

impl ValueEnum for SimilarityMethod {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            SimilarityMethod::Kmer,
            SimilarityMethod::Minimizer,
            SimilarityMethod::Strobemer,
        ]
    }

    fn to_possible_value<'a>(&self) -> Option<PossibleValue> {
        Some(match self {
            SimilarityMethod::Kmer => PossibleValue::new("kmer"),
            SimilarityMethod::Minimizer => PossibleValue::new("minimizer"),
            SimilarityMethod::Strobemer => PossibleValue::new("strobemer"),
        })
    }
}

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
        let estimated_distance: f64 = match args.method {
            SimilarityMethod::Kmer => kmer_similarity(
                &base_seq,
                &mod_seq,
                args.distance_function.clone(),
                args.k.clone(),
                args.step.clone(),
            )?,
            SimilarityMethod::Minimizer => minimizer_similarity(
                &base_seq,
                &mod_seq,
                args.distance_function.clone(),
                args.k.clone(),
                args.minimizer_window_length.clone(),
                args.step.clone(),
            )?,
            SimilarityMethod::Strobemer => strobemer_similarity(
                &base_seq,
                &mod_seq,
                args.distance_function.clone(),
                args.strobemer_order.clone(),
                args.strobemer_length.clone(),
                args.strobemer_window_gap.clone(),
                args.strobemer_window_length.clone(),
                args.step.clone(),
            )?,
        };

        edit_distance_sums
            .entry(record.edit_distance.clone())
            .and_modify(|v| *v += estimated_distance)
            .or_insert(0.0);
        edit_distance_squared_sums
            .entry(record.edit_distance.clone())
            .and_modify(|v| *v += estimated_distance * estimated_distance)
            .or_insert(0.0);
        edit_distance_counts
            .entry(record.edit_distance.clone())
            .and_modify(|v| *v += 1.)
            .or_insert(0.0);

        print!("\rString pair #{:?} has been processed!", i);
        io::stdout().flush()?;
    }

    // Compute mean and confidence intervals from rolling data.
    let mut lower_bounds: HashMap<usize, f64> = HashMap::new();
    let mut upper_bounds: HashMap<usize, f64> = HashMap::new();
    let mut means: HashMap<usize, f64> = HashMap::new();
    for (edit_distance, count) in &edit_distance_counts {
        let sum = edit_distance_sums.get(&edit_distance).unwrap();
        let squared_sum =
            edit_distance_squared_sums.get(&edit_distance).unwrap();

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

    let mut sorted_edit_distances: Vec<&usize> =
        edit_distance_counts.keys().collect();
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
