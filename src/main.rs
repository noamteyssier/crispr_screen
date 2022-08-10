use std::path::Path;
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
use utils::io::load_dataframe;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {

    /// Filepath of the input count matrix
    #[clap(short, long, value_parser)]
    input: String,

    /// Labels for Control Samples
    #[clap(short, long, value_parser, required=true)]
    controls: Vec<String>,

    /// Labels for Treatment Samples
    #[clap(short, long, value_parser, required=true)]
    treatments: Vec<String>,

    /// Output Prefix
    #[clap(short, long, value_parser, default_value="results")]
    output: String,

    /// Normalization Option
    #[clap(short, long, value_parser, default_value="median")]
    norm: String,

    /// Aggregation Option
    #[clap(short='g', long, value_parser, default_value="rra")]
    agg: String,

    /// Permutations
    #[clap(short, long, value_parser, default_value="100")]
    permutations: usize,

    /// Alpha Threshold
    #[clap(short, long, value_parser, default_value="0.1")]
    alpha: f64,

    /// Non-Targeting Control Token
    #[clap(short, long, value_parser, default_value="non-targeting")]
    ntc_token: String
}

fn main() {
    let args = Args::parse();

    // validate input path
    let path = if Path::new(&args.input).exists() { 
        args.input 
    } else { 
        panic!("Provided Input Does Not Exist: {}", args.input) 
    };

    // validate normalization method
    let norm_method = match args.norm.as_str() {
        "median" => Normalization::MedianRatio,
        "total" => Normalization::Total,
        _ => panic!("Unexpected normalization method provided: {}", args.norm)
    };

    // validate aggregation method
    let agg = match args.agg.as_str() {
        "rra" => GeneAggregation::AlpaRRA { alpha: args.alpha, npermutations: args.permutations },
        "inc" => GeneAggregation::Inc { token: &args.ntc_token },
        _ => panic!("Unexpected aggregation method provided: {}", args.agg)
    };

    let labels_controls = args.controls;
    let labels_treatments = args.treatments;
    let frame = load_dataframe(&path).unwrap();

    mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        &args.output,
        &norm_method,
        &agg
    ).unwrap();
}
