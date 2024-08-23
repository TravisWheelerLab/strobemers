use std::fs::File;
use std::path::Path;
use std::io::{BufReader, Write};

use anyhow::Result;
use clap::Parser;
use bio::io::fasta::Reader;


#[derive(Debug, Parser)]
#[command(about, author, version)]
pub struct CommonArgs {
    #[arg(required = true, value_name = "REF FILE", help = "References FASTA file")]
    pub references_file: String,
    #[arg(required = true, value_name = "QUERY FILE", help = "Query FAST(A/Q) file")]
    pub query_file: String,
    //#[arg(short='m', value_name = "METHOD()", help = "E.g. jaccard_similarity")] // unused
    //pub similarity_method: String,

    #[arg(short='e', value_name = "STRING", help="experiment name: output subdirectory to put CSV file")]
    pub experiment_name: String
}

fn create_csv(args_common: &CommonArgs, file_name: &String) -> Result<File>{
    let file_path = Path::new(&std::env::var("CARGO_MANIFEST_DIR")?)
        .join("data")
        .join(args_common.experiment_name.clone());
    if !std::path::Path::new(&file_path).exists() {
        std::fs::create_dir(&file_path)?;
    }
    let file = File::create(&file_path.join(file_name))?;
    Ok(file)
}

pub fn create_estimation_csv_with_headers(args_common: &CommonArgs) -> Result<File>{
    let mut file = create_csv(args_common, &"estimation-output.csv".to_string())?;
    writeln!(file, "ref_name,query_name,seed_name,estimation,ref_time,query_time")?;
    Ok(file)
}

pub fn create_alignment_csv_with_headers(args_common: &CommonArgs) -> Result<File>{
    let mut file = create_csv(args_common, &"alignment-output.csv".to_string())?;
    writeln!(file, "ref_name,query_name,seed_name,edit_distance,edit_distance_time")?;
    Ok(file)
}

pub fn create_query_reader(args_common: &CommonArgs) -> Result<Reader<BufReader<File>>> {
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = Reader::from_file(
        Path::new(&project_dir)
            .join("data")
            .join(&args_common.experiment_name)
            .join(&args_common.query_file)
        )?;
    Ok(query_reader)
}

pub fn create_reference_reader(args_common: &CommonArgs) -> Result<Reader<BufReader<File>>> {
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = Reader::from_file(
        Path::new(&project_dir)
            .join("data")
            .join(&args_common.experiment_name)
            .join(&args_common.references_file)
        )?;
    Ok(query_reader)
}

pub fn print_update(i: usize) {
    print!("\r comparisons done: {:?}", i);
    std::io::stdout().flush().unwrap();
}