use anyhow::Result;
use clap::Parser;

use alignment_free_methods::cli::*;

// --------------------------------------------------

#[derive(Debug, Parser)]
pub struct StrobemerArgs {
    #[command(flatten)]
    pub common_args: CommonArgs,    
    #[command(flatten)]
    pub strobemer_args: StrobemerSpecificArgs
}

// --------------------------------------------------

fn main() -> Result<()>{
    let args = StrobemerArgs::parse();
    match run(&args.common_args, &args.strobemer_args) {
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
        Ok(results_to_save) => {
            let seed_name = args.strobemer_args.repr();
            args.common_args.save_results_to_csv(&seed_name, &results_to_save)?;
        }
    }
    Ok(())
}