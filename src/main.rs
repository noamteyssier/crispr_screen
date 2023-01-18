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
use aggregation::GeneAggregation;
use utils::{io::load_dataframe, logging::Logger};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {

    /// Filepath of the input count matrix
    #[arg(short, long)]
    input: String,

    /// Labels for Control Samples
    #[arg(short, long, use_value_delimiter=true, value_delimiter = ' ')]
    controls: Vec<String>,

    /// Labels for Treatment Samples
    #[arg(short, long, use_value_delimiter=true, value_delimiter = ' ')]
    treatments: Vec<String>,

    /// Output Prefix
    #[arg(short, long, default_value="results")]
    output: String,

    /// Normalization Option
    #[arg(short, long, default_value="median-ratio")]
    // norm: String,
    norm: Normalization,

    /// Aggregation Option
    #[arg(short='g', long, default_value="rra")]
    agg: String,

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

    /// Multiple Hypothesis Correction (bonferroni, bh, by)
    #[arg(short='f', long, default_value="bh")]
    correction: String
}

fn main() -> Result<()> {
    let args = Args::parse();

    // validate input path
    let path = if Path::new(&args.input).exists() { 
        args.input 
    } else { 
        panic!("Provided Input Does Not Exist: {}", args.input) 
    };

    // validate aggregation method
    let agg = match args.agg.as_str() {
        "rra" => GeneAggregation::AlpaRRA { alpha: args.alpha, npermutations: args.permutations },
        "inc" => GeneAggregation::Inc { token: &args.ntc_token },
        _ => panic!("Unexpected aggregation method provided: {}", args.agg)
    };

    // create logger based on quiet option
    let logger = if args.quiet {
        Logger::new_silent()
    } else {
        Logger::new()
    };

    // create multiple hypothesis correction from option
    let correction = match args.correction.as_str() {
        "bonferroni" => Procedure::Bonferroni,
        "bh" | "fdr" => Procedure::BenjaminiHochberg,
        "by" => Procedure::BenjaminiYekutieli,
        _ => panic!("Unexpected correction method provided: {}", args.correction)
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

