use std::fs::File;
use std::path::Path;
use std::io::{BufReader, Write};
use std::time::Instant;

use anyhow::Result;
use clap::Parser;
use bio::io::fasta::Reader;

use crate::*;

// --------------------------------------------------
// Arguments.

// Common arguments are positional arguments for input and output file names.
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

#[derive(Debug, Parser)]
pub struct KmerSpecificArgs {
    #[arg(short='k', value_name = "INT", help="Substring length")]
    pub k: usize,

    #[arg(long="ref-index", value_name = "REF", help="w-max: window selection parameter", default_value_t=1)]
    pub ref_index: usize,
}

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
// Helper functions for run() -- mostly code deduplication.

pub fn create_reader(file_path_from_manifest: &str) -> Result<Reader<BufReader<File>>>{
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = Reader::from_file(
        Path::new(&project_dir).join(file_path_from_manifest)
    )?;
    Ok(query_reader)
}

// Does what it says. Uses carriage return.
pub fn print_update(i: usize) {
    print!("\r comparisons done: {:?}", i);
    std::io::stdout().flush().unwrap();
}

pub fn generic_seed_comparison<T: SeedTraits>(args: &CommonArgs, seed_specific_args: &T) -> Result<Vec<String>> {
    let mut results_to_save = Vec::new();

    let query_reader = create_reader(&args.query_file)?;
    for query_record in query_reader.records() {
        // There should only be 1 query record!
        let query_record = query_record?;
        let query_seed_start_time = Instant::now();
        
        let query_seeds = seed_specific_args.generate_seeds(query_record.seq())?;
        let query_time = query_seed_start_time.elapsed().as_secs_f64();

        let references_reader = create_reader(&args.references_file)?;
        let mut i: usize = 0;
        for reference_record in references_reader.records() {
            let reference_record = reference_record?;
            let reference_time = Instant::now();
            let ref_seeds = seed_specific_args.generate_seeds(reference_record.seq())?;
            let reference_time = reference_time.elapsed().as_secs_f64();
            let estimation = jaccard_similarity(&ref_seeds,&query_seeds)?;
            print_update(i);
            i += 1;

            let result = format!("{},{},{}",
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
// SeedTraits allows run() to work with arbitrary [seed]SpecificArgs.
// All you need to do is define how generate_seeds works for that kind of seed.
// There is also a repr() function which helps with more code deduplication.
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

// --------------------------------------------------
// More deduplication -- SaveToCVS makes it super easy to save results.
// You need a results_to_save object and seed_name though. results_to_save is just a Vec<String>
// containing "{estimation},{query_time},{ref_time}" lines. This trait takes care of the rest.
pub trait SaveToCSV {
    fn save_results_to_csv(&self, seed_name: &String, results_to_save: &Vec<String>) -> Result<()>;
}

impl SaveToCSV for CommonArgs{
    fn save_results_to_csv(&self, seed_name: &String, results_to_save: &Vec<String>) -> Result<()>{
        let mut file = File::create(&self.output_file)?;
        writeln!(file, "seed_name,query_name,ref_name,estimation,query_time,ref_time")?;
        for result in results_to_save {
            writeln!(file, "{},{},{},{}",
                seed_name,
                &self.query_file,
                &self.references_file,
                result
            )?;
        }    
        Ok(())
    }        
}