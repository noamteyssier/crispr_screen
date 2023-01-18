use anyhow::Result;
use std::path::Path;
use adjustp::Procedure;
use clap::Parser;

mod aggregation;
mod model;
mod enrich;
mod norm;
mod utils;
mod differential_expression;

use differential_expression::mageck;
use norm::Normalization;
use aggregation::{GeneAggregation, GeneAggregationSelection};
use utils::{io::load_dataframe, logging::Logger, Adjustment};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {

    /// Filepath of the input count matrix
    #[arg(short, long)]
    input: String,

    /// Labels for Control Samples
    #[arg(short, long, num_args=1.., required=true)]
    controls: Vec<String>,

    /// Labels for Treatment Samples
    #[arg(short, long, num_args=1.., required=true)]
    treatments: Vec<String>,

    /// Output Prefix
    #[arg(short, long, default_value="results")]
    output: String,

    /// Normalization Option
    #[arg(short, long, default_value="median-ratio")]
    norm: Normalization,

    /// Aggregation Option
    #[arg(short='g', long, default_value="rra")]
    agg: GeneAggregationSelection,

    /// Permutations
    #[arg(short, long, default_value="100")]
    permutations: usize,

    /// Alpha Threshold
    #[arg(short, long, default_value="0.1")]
    alpha: f64,

    /// Non-Targeting Control Token
    #[arg(long, default_value="non-targeting")]
    ntc_token: String,

    /// Do not write logging information
    #[arg(short, long)]
    quiet: bool,

    /// Multiple Hypothesis Correction
    #[arg(short='f', long, default_value="bh")]
    correction: Adjustment
}

fn main() -> Result<()> {
    let args = Args::parse();

    // validate input path
    let path = if Path::new(&args.input).exists() { 
        args.input 
    } else { 
        panic!("Provided Input Does Not Exist: {}", args.input) 
    };

    // assign and parameterize gene aggregation method
    let agg = match args.agg {
        GeneAggregationSelection::RRA => {
            GeneAggregation::AlpaRRA { 
                alpha: args.alpha, 
                npermutations: args.permutations, 
            }
        },
        GeneAggregationSelection::Inc => {
            GeneAggregation::Inc {
                token: &args.ntc_token,
            }
        }
    };

    // create logger based on quiet option
    let logger = if args.quiet {
        Logger::new_silent()
    } else {
        Logger::new()
    };

    // create multiple hypothesis correction from option
    let correction = match args.correction {
        Adjustment::Bf => Procedure::Bonferroni,
        Adjustment::Bh => Procedure::BenjaminiHochberg,
        Adjustment::By => Procedure::BenjaminiYekutieli
    };

    let labels_controls = args.controls;
    let labels_treatments = args.treatments;
    let frame = load_dataframe(&path).unwrap();

    let mageck_results = mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        &args.output,
        &args.norm,
        &agg,
        &logger,
        &correction
    );

    match mageck_results {
        Err(e) => {
            println!("ERROR: {}", e);
            Ok(())
        },
        Ok(_) => {
            Ok(())
        }
    }
}

