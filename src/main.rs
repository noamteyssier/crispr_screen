use adjustp::Procedure;
use anyhow::Result;
use clap::Parser;
use cli::Cli;
use io::SimpleFrame;
use std::path::Path;

pub mod aggregation;
pub mod cli;
pub mod differential_expression;
pub mod enrich;
pub mod io;
pub mod model;
pub mod norm;
pub mod utils;

use aggregation::{GeneAggregation, GeneAggregationSelection};
use differential_expression::mageck;
use utils::{config::Configuration, logging::Logger, Adjustment};

fn main() -> Result<()> {
    let args = Cli::parse();

    // validate input path
    let path = if Path::new(&args.input).exists() {
        args.input
    } else {
        panic!("Provided Input Does Not Exist: {}", args.input)
    };

    // set rayon threads
    if let Some(t) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
    }

    // assign and parameterize gene aggregation method
    let agg = match args.agg {
        GeneAggregationSelection::RRA => GeneAggregation::AlpaRRA {
            alpha: args.alpha,
            npermutations: args.permutations,
            adjust_alpha: !args.no_adjust_alpha,
            fdr: args.fdr,
        },
        GeneAggregationSelection::Inc => GeneAggregation::Inc {
            token: &args.ntc_token,
            fdr: args.fdr,
            group_size: args.inc_group_size,
            use_product: args.inc_product,
            n_draws: args.n_draws,
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
        args.seed,
        &args.output,
    );

    let labels_controls = args.controls;
    let labels_treatments = args.treatments;
    let frame = SimpleFrame::from_filepath(&path)?;

    let mageck_results = mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        &config,
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
