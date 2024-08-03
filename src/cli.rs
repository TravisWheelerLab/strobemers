use std::path::PathBuf;
use anyhow::Result;
use clap::Parser;
use rusqlite::{self, Connection};


#[derive(Debug, Parser)]
#[command(about, author, version)]
pub struct CommonArgs {
    #[arg(required = true, value_name = "REF FILE", help = "References FASTA file")]
    pub references_file: String,
    #[arg(required = true, value_name = "QUERY FILE", help = "Query FAST(A/Q) file")]
    pub query_file: String,
    //#[arg(short='m', value_name = "METHOD()", help = "E.g. jaccard_similarity")] // unused
    //pub similarity_method: String,
    #[arg(short='s', default_value_t = 1, value_name = "INT PARAM")]
    pub step: usize
}


pub fn initialize_comparison_db(filename: PathBuf) -> Result<Connection> {
    let conn = Connection::open(filename)?;
    conn.execute(
        "CREATE TABLE IF NOT EXISTS comparisons (
            query_name TEXT NOT NULL,
            reference_name TEXT NOT NULL,
            seed_name TEXT NOT NULL,
            score REAL NOT NULL,
            time REAL NOT NULL,
            UNIQUE(query_name, reference_name, seed_name)
        )",
        [],
    )?;
    Ok(conn)
}