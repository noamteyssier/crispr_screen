use adjustp::Procedure;
use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands, DiffAbundanceArgs, IncArgs, InputArgs, MiscArgs, RraArgs, SgrnaColumns};
use io::SimpleFrame;
use run_aggregation::run_aggregation;
use std::path::Path;

pub mod aggregation;
pub mod cli;
pub mod differential_expression;
pub mod enrich;
pub mod io;
pub mod model;
pub mod norm;
pub mod run_aggregation;
pub mod utils;

use aggregation::{GeneAggregation, GeneAggregationSelection};
use differential_expression::mageck;
use utils::{config::Configuration, logging::Logger, Adjustment};

#[allow(clippy::too_many_arguments)]
fn test(
    input_args: InputArgs,
    prefix: String,
    diff_args: DiffAbundanceArgs,
    agg: GeneAggregationSelection,
    rra: RraArgs,
    inc: IncArgs,
    misc: MiscArgs,
    skip_agg: bool,
) -> Result<()> {
    // validate input path
    let path = if Path::new(&input_args.input).exists() {
        input_args.input
    } else {
        panic!("Provided Input Does Not Exist: {}", input_args.input)
    };

    // set rayon threads
    if let Some(t) = misc.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
    }

    // assign and parameterize gene aggregation method
    let agg = match agg {
        GeneAggregationSelection::RRA => GeneAggregation::AlpaRRA {
            alpha: rra.alpha,
            npermutations: rra.permutations,
            adjust_alpha: !rra.no_adjust_alpha,
            fdr: misc.fdr,
        },
        GeneAggregationSelection::Inc => GeneAggregation::Inc {
            token: &inc.ntc_token,
            group_size: inc.inc_group_size,
            use_product: inc.inc_product,
            n_draws: inc.n_draws,
            fdr: misc.fdr,
        },
    };

    // create logger based on quiet option
    let logger = if misc.quiet {
        Logger::new_silent()
    } else {
        Logger::new()
    };

    // create multiple hypothesis correction from option
    let correction = match misc.correction {
        Adjustment::Bf => Procedure::Bonferroni,
        Adjustment::Bh => Procedure::BenjaminiHochberg,
        Adjustment::By => Procedure::BenjaminiYekutieli,
    };

    let config = Configuration::new(
        diff_args.norm,
        agg,
        correction,
        diff_args.model_choice,
        misc.seed,
        &prefix,
    );

    let labels_controls = input_args.controls;
    let labels_treatments = input_args.treatments;
    let frame = SimpleFrame::from_filepath(&path)?;

    let mageck_results = mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        &config,
        &logger,
        skip_agg,
    );

    match mageck_results {
        Err(e) => {
            println!("ERROR: {e}");
            Ok(())
        }
        Ok(_) => Ok(()),
    }
}

fn aggregate(
    input: String,
    prefix: String,
    columns: SgrnaColumns,
    agg: GeneAggregationSelection,
    rra: RraArgs,
    inc: IncArgs,
    misc: MiscArgs,
) -> Result<()> {
    // validate input path
    let path = if Path::new(&input).exists() {
        input
    } else {
        panic!("Provided Input Does Not Exist: {}", input)
    };

    // set rayon threads
    if let Some(t) = misc.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
    }

    // assign and parameterize gene aggregation method
    let agg = match agg {
        GeneAggregationSelection::RRA => GeneAggregation::AlpaRRA {
            alpha: rra.alpha,
            npermutations: rra.permutations,
            adjust_alpha: !rra.no_adjust_alpha,
            fdr: misc.fdr,
        },
        GeneAggregationSelection::Inc => GeneAggregation::Inc {
            token: &inc.ntc_token,
            group_size: inc.inc_group_size,
            use_product: inc.inc_product,
            n_draws: inc.n_draws,
            fdr: misc.fdr,
        },
    };

    // create logger based on quiet option
    let logger = if misc.quiet {
        Logger::new_silent()
    } else {
        Logger::new()
    };

    // create multiple hypothesis correction from option
    let correction = match misc.correction {
        Adjustment::Bf => Procedure::Bonferroni,
        Adjustment::Bh => Procedure::BenjaminiHochberg,
        Adjustment::By => Procedure::BenjaminiYekutieli,
    };

    let config = Configuration::new_agg(agg, correction, misc.seed, &prefix);
    let frame = SimpleFrame::from_filepath(&path)?;

    run_aggregation(&frame, columns, &config, &logger)
}

fn main() -> Result<()> {
    let args = Cli::parse();
    match args.command {
        Commands::Test {
            input,
            prefix,
            diff_args,
            agg,
            rra,
            inc,
            misc,
            skip_agg,
        } => test(input, prefix, diff_args, agg, rra, inc, misc, skip_agg),
        Commands::Agg {
            input,
            prefix,
            columns,
            agg,
            rra,
            inc,
            misc,
        } => aggregate(input, prefix, columns, agg, rra, inc, misc),
    }
}
