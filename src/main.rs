use adjustp::Procedure;
use anyhow::Result;
use clap::Parser;
use io::SimpleFrame;
use model::ModelChoice;
use std::path::Path;

mod aggregation;
mod differential_expression;
mod enrich;
mod io;
mod model;
mod norm;
mod utils;

use aggregation::{GeneAggregation, GeneAggregationSelection};
use differential_expression::mageck;
use norm::Normalization;
use utils::{logging::Logger, Adjustment, config::Configuration};

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
    #[arg(short, long, default_value = "results")]
    output: String,

    /// Normalization Option
    #[arg(short, long, default_value = "median-ratio")]
    norm: Normalization,

    /// Aggregation Option
    #[arg(short = 'g', long, default_value = "rra")]
    agg: GeneAggregationSelection,

    /// Permutations
    #[arg(short, long, default_value = "100")]
    permutations: usize,

    /// Alpha Threshold
    #[arg(short, long, default_value = "0.25")]
    alpha: f64,

    /// Do not adjust alpha threshold for RRA.
    #[arg(long)]
    no_adjust_alpha: bool,

    /// Non-Targeting Control Token
    #[arg(long, default_value = "non-targeting")]
    ntc_token: String,

    /// Do not write logging information
    #[arg(short, long)]
    quiet: bool,

    /// Multiple Hypothesis Correction
    #[arg(short = 'f', long, default_value = "bh")]
    correction: Adjustment,

    /// Least Squares Model Choice
    #[arg(short, long, default_value = "wols")]
    model_choice: ModelChoice,
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
        GeneAggregationSelection::RRA => GeneAggregation::AlpaRRA {
            alpha: args.alpha,
            npermutations: args.permutations,
            adjust_alpha: !args.no_adjust_alpha,
        },
        GeneAggregationSelection::Inc => GeneAggregation::Inc {
            token: &args.ntc_token,
        },
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
        Adjustment::By => Procedure::BenjaminiYekutieli,
    };

    let config = Configuration::new(
        args.norm, 
        agg, 
        correction, 
        args.model_choice,
        &args.output,
    );

    let labels_controls = args.controls;
    let labels_treatments = args.treatments;
    let frame = SimpleFrame::from_filepath(&path)?;

    let mageck_results = mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        config,
        &logger,
    );

    match mageck_results {
        Err(e) => {
            println!("ERROR: {e}");
            Ok(())
        }
        Ok(_) => Ok(()),
    }
}
