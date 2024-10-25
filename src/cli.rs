//! Command-line interface code.

use std::fs::File;
use std::io::{BufReader, BufRead, Write};
use std::time::Instant;

use anyhow::Result;
use clap::Parser;

use crate::*;

// --------------------------------------------------
// Arguments.

/// Positional arguments for input and output file names.
#[derive(Debug, Parser)]
#[command(about, author, version)]
pub struct CommonArgs {
    // Everything is relative to the root of the directory!
    #[arg(required = true, value_name = "REF", help = "References .FASTA file")]
    pub references_file: String,
    #[arg(required = true, value_name = "QUERY", help = "Query .FAST(A/Q) file")]
    pub query_file: String,
    #[arg(required = true, value_name = "OUT", help = "Output .csv file")]
    pub output_file: String,
}

/// Keyword arguments for k-mers.
#[derive(Debug, Parser)]
pub struct KmerSpecificArgs {
    #[arg(short='k', value_name = "INT", help="Substring length")]
    pub k: usize,

    #[arg(long="ref-index", value_name = "REF", help="w-max: window selection parameter", default_value_t=1)]
    pub ref_index: usize,
}

/// Keyword arguments for Strobemers
#[derive(Debug, Parser)]
pub struct StrobemerSpecificArgs {
    #[arg(short='p', value_name = "PROTOCOL", help="Strobemer generation protocol", default_value = "rand")]
    pub protocol: Protocol,
    #[arg(short='o', value_name = "ORDER", help="Strobemer order/how many strobes", default_value_t = 2)]
    pub order: usize,
    #[arg(short='l', value_name = "L", help="Strobe length")]
    pub strobe_length: usize,
    #[arg(long="w-min", value_name = "W_MIN", help="w-min: window selection parameter")]
    pub w_min: usize,
    #[arg(long="w-max", value_name = "W_MAX", help="w-max: window selection parameter")]
    pub w_max: usize,

    #[arg(long="ref-index", value_name = "REF", help="w-max: window selection parameter", default_value_t=1)]
    pub ref_index: usize,
}

// --------------------------------------------------
// Helper functions

fn parse_fasta(filename: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut current_sequence = Vec::new();
    let mut current_id = String::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            // New sequence started, store previous sequence if any
            if !current_sequence.is_empty() {
                sequences.push((current_id, current_sequence.clone()));
                current_sequence.clear();
            }
            // id to new sequence
            current_id = line[1..].to_string();
        } else {
            // Process sequence data
            let converted_sequence: Vec<u8> = line
                .chars()
                .map(|c| match c {
                    'A' => 0,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    'a' => 0,
                    'c' => 1,
                    'g' => 2,
                    't' => 3,
                    _ => 0, // For all other IUB codes
                })
                .collect();
            current_sequence.extend(converted_sequence);
        }
    }

    // Store the last sequence
    if !current_sequence.is_empty() {
        sequences.push((current_id, current_sequence));
    }

    Ok(sequences)
}

/// Prints an update to stdout with carriage return.
pub fn print_update(i: usize) {
    print!("\r comparisons done: {:?}", i);
    std::io::stdout().flush().unwrap();
}

/// Generates a [`Vec`]\<[`SeedObject`]\> based on how seed_specific_args implements
/// [`SeedTraits`] (i.e. [`StrobemerSpecificArgs`] will generate Strobemers).
pub fn generic_seed_comparison<T: SeedTraits>(args: &CommonArgs, seed_specific_args: &T) -> Result<Vec<String>> {
    let mut results_to_save = Vec::new();

    for (query_id, query_seq) in parse_fasta(&args.query_file)? {
        let query_seed_start_time = Instant::now();
        
        let query_seeds = seed_specific_args.generate_seeds(&query_seq)?;
        let query_time = query_seed_start_time.elapsed().as_secs_f64();

        let mut i: usize = 0;
        for (reference_id, reference_seq) in parse_fasta(&args.references_file)? {
            let reference_time = Instant::now();
            let ref_seeds = seed_specific_args.generate_seeds(&reference_seq)?;
            let reference_time = reference_time.elapsed().as_secs_f64();
            let estimation = jaccard_similarity(&ref_seeds,&query_seeds)?;
            print_update(i);
            i += 1;

            let result = format!("{},{},{},{},{},{}",
                seed_specific_args.repr(),
                query_id,
                reference_id,
                estimation,
                query_time,
                reference_time
            );
            results_to_save.push(result);
        }
    }
    Ok(results_to_save)
}

// --------------------------------------------------
/// Allows code to take a generic type and generate a [`Vec`]\<[`SeedObject`].
pub trait SeedTraits {
    fn generate_seeds(&self, seq: &[u8]) -> Result<Vec<SeedObject>>;
    fn repr(&self) -> String;
}

impl SeedTraits for KmerSpecificArgs {
    fn generate_seeds(&self, seq: &[u8]) -> Result<Vec<SeedObject>>{
        seq_to_kmers(seq, &self)
    }
    fn repr(&self) -> String {
        format!("{}-mers", {self.k})
    }
}

impl SeedTraits for StrobemerSpecificArgs {
    fn generate_seeds(&self, seq: &[u8]) -> Result<Vec<SeedObject>>{
        seq_to_strobemers(seq, &self)
    }
    fn repr(&self) -> String {        
        format!("({}.{}.{}.{})-{}strobemers",
            self.order,
            self.strobe_length,
            self.w_min,
            self.w_max,
            self.protocol
        )
    }
}


/// A trait to make saving results easier.
pub trait SaveToCSV {
    fn save_results_to_csv(&self, results_to_save: &Vec<String>) -> Result<()>;
}

impl SaveToCSV for CommonArgs{
    fn save_results_to_csv(&self, results_to_save: &Vec<String>) -> Result<()>{
        let mut file = File::create(&self.output_file)?;
        writeln!(file, "seed_name,query_name,ref_name,estimation,query_time,ref_time")?;
        for result in results_to_save {
            writeln!(file, "{}", result)?;
        }    
        Ok(())
    }        
}