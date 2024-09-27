use std::fs::File;
use std::path::Path;
use std::io::{BufReader, Write};
use std::str::FromStr;
use std::fmt;

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
}

#[derive(Debug, Clone)]
pub enum Protocol {
    Rand,
    Min,
    Hybrid,
}

impl FromStr for Protocol {
    type Err = String;

    fn from_str(input: &str) -> Result<Protocol, Self::Err> {
        match input {
            "rand" => Ok(Protocol::Rand),
            "min" => Ok(Protocol::Min),
            "hybrid" => Ok(Protocol::Hybrid),
            _ => Err(format!("Invalid protocol: {}", input)),
        }
    }
}

impl fmt::Display for Protocol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let protocol_str = match self {
            Protocol::Rand => "rand",
            Protocol::Min => "min",
            Protocol::Hybrid => "hybrid",
        };
        write!(f, "{}", protocol_str)
    }
}

#[derive(Debug, Parser)]
pub struct StrobemerArgs {
    #[command(flatten)]
    pub common: CommonArgs,

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

pub fn create_reader(file_path_from_manifest: &str) -> Result<Reader<BufReader<File>>>{
    let project_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    let query_reader = Reader::from_file(
        Path::new(&project_dir).join(file_path_from_manifest)
    )?;
    Ok(query_reader)
}

pub fn print_update(i: usize) {
    print!("\r comparisons done: {:?}", i);
    std::io::stdout().flush().unwrap();
}

pub fn format_strobemer_seed_name(args: &StrobemerArgs) -> String {
    let seed_name = format!("({}.{}.{}.{})-{}strobemers",
        &args.order,
        &args.strobe_length,
        &args.w_min,
        &args.w_max,
        &args.protocol
    );
    seed_name
}