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
    #[clap(short, long, value_parser)]
    norm: Option<String>,

    /// Permutations
    #[clap(short, long, value_parser, default_value="100")]
    permutations: usize,

    /// Alpha Threshold
    #[clap(short, long, value_parser, default_value="0.1")]
    alpha: f64
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
    let norm_method = match args.norm {
        Some(x) => match x.as_str() {
            "median" => Normalization::MedianRatio,
            "total" => Normalization::Total,
            _ => panic!("Unexpected normalization method provided: {}", x)
        },
        None => Normalization::MedianRatio,
    };
    let labels_controls = args.controls;
    let labels_treatments = args.treatments;
    let frame = load_dataframe(&path).unwrap();
    let agg = GeneAggregation::AlpaRRA { alpha: args.alpha, npermutations: args.permutations };
    mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        &args.output,
        norm_method,
        agg
    ).unwrap();
}
