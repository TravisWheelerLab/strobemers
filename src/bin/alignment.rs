use std::time::Instant;
use bio::alignment::distance::levenshtein;
use anyhow::Result;
use clap::Parser;
use glob::glob;
use std::path::Path;

use alignment_free_methods::cli::*;

// --------------------------------------------------
fn main() {
    if let Err(e) = run(CommonArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
// See this repo's README file for pseudocode
fn run(args: CommonArgs) -> Result<()> {
    let project_dir = Path::new(&std::env::var("CARGO_MANIFEST_DIR")?);

    let query_files = glob(
        &format!("{:?}/data/experiment{}/{}", project_dir, &args.experiment_name, "query*.csv")
    )?;

    for query_filename in query_files {
        let mut csv_file = create_csv_with_alignment_headers(&args, )?;
        let seed_name = format!("alignment");
        let query_reader = create_query_reader(&args)?;

        let mut i = 0;
        for query_record in query_reader.records() {
            let query_record = query_record?;
            let reference_reader = create_reference_reader(&args)?;
            for reference_record in reference_reader.records() {
                i += 1;
                let reference_record = reference_record?;

                print_update(i);

                let start = Instant::now();
                let estimation = levenshtein(query_record.seq(),reference_record.seq());
                let duration = start.elapsed().as_secs_f64();

                write_to_alignment_csv(&mut csv_file,
                    &reference_record.id(),
                    &query_record.id(),
                    &seed_name,
                    &estimation.to_string(),
                    &duration.to_string()
                )?
            
            }
        }
    }
    Ok(())

}