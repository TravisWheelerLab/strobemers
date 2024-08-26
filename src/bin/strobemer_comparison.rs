use std::time::Instant;
use anyhow::Result;
use bio::alignment::distance::levenshtein;
use clap::Parser;

use alignment_free_methods::cli::*;
use alignment_free_methods::{seq_to_strobemers, jaccard_similarity};

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
    #[arg(long="w-min", value_name = "INT", help="w-min: window selection parameter")]
    w_min: usize,
    #[arg(long="w-max", value_name = "INT", help="w-max: window selection parameter")]
    w_max: usize,
    #[arg(long="run-alignment", default_value_t = false)]
    run_alignment: bool,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(StrobemerArgs::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run(args: StrobemerArgs) -> Result<()> {
    let mut csv_file = match args.run_alignment {
        true => create_csv_with_estimation_and_alignment_headers(&args.common)?,
        false => create_csv_with_estimation_headers(&args.common)?
    };

    let seed_name = format!("({}.{}.{}.{})-{}strobemers",
        &args.order,
        &args.strobe_length,
        &args.w_min,
        &args.w_max,
        &args.protocol
    );

    let query_reader = create_query_reader(&args.common)?;

    let mut i = 0;
    for query_record in query_reader.records() {
        let query_record = query_record?;
        let query_seed_start_time = Instant::now();
        let query_seeds = seq_to_strobemers(
            query_record.seq(),
            args.order.clone(),
            args.strobe_length.clone(),
            args.w_min.clone(),
            args.w_max.clone(),
            1,
            &args.protocol
        )?;
        let query_time = query_seed_start_time.elapsed().as_secs_f64();

        let reference_reader = create_reference_reader(&args.common)?;
        for reference_record in reference_reader.records() {
            i += 1;
            let reference_record = reference_record?;

            let reference_time = Instant::now();
            let reference_seeds = seq_to_strobemers(
                reference_record.seq(),
                args.order.clone(),
                args.strobe_length.clone(),
                args.w_min.clone(),
                args.w_max.clone(),
                1,
                &args.protocol
            )?;
            let reference_time = reference_time.elapsed().as_secs_f64();

            let estimation = jaccard_similarity(&reference_seeds,&query_seeds,)?;


            print_update(i);

            match args.run_alignment {
                true => {
                    let alignment_time = Instant::now();
                    let edit_distance = levenshtein(reference_record.seq(), query_record.seq());
                    let alignment_time = alignment_time.elapsed().as_secs_f64();
                    write_to_estimation_and_alignment_csv(&mut csv_file,
                        reference_record.id(),
                        query_record.id(),
                        &seed_name,
                        &estimation.to_string(),
                        &reference_time.to_string(),
                        &query_time.to_string(),
                        &edit_distance.to_string(),
                        &alignment_time.to_string(),
                    )?;
                },
                false => write_to_estimation_csv(&mut csv_file,
                    reference_record.id(),
                    query_record.id(),
                    &seed_name,
                    &estimation.to_string(),
                    &reference_time.to_string(),
                    &query_time.to_string()
                )?
            };
            

            
        }
    }
    Ok(())
}