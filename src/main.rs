use adjustp::Procedure;
use anyhow::Result;
use clap::Parser;
use cli::{
    Cli, Commands, DiffAbundanceArgs, GeopaggArgs, IncArgs, InputArgs, MiscArgs, RraArgs,
    SgrnaColumns,
};
use geopagg::WeightConfig;
use regex::Regex;
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

use aggregation::{GeneAggregation, GeneAggregationSelection, GeoPAGGWeightConfigEnum};
use differential_expression::mageck;
use io::load_dataframe;
use utils::{config::Configuration, logging::Logger, Adjustment};

#[allow(clippy::too_many_arguments)]
fn test(
    input_args: InputArgs,
    prefix: String,
    diff_args: DiffAbundanceArgs,
    agg: GeneAggregationSelection,
    rra: RraArgs,
    inc: IncArgs,
    geopagg: GeopaggArgs,
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
            token: &misc.ntc_token,
            group_size: inc.inc_group_size,
            use_product: inc.inc_product,
            n_draws: inc.n_draws,
            fdr: misc.fdr,
        },
        GeneAggregationSelection::GeoPAGG => GeneAggregation::GeoPAGG {
            token: &misc.ntc_token,
            fdr: misc.fdr,
            weight_config: {
                match geopagg.weight_config {
                    GeoPAGGWeightConfigEnum::DropFirst => WeightConfig::DropFirst {
                        alpha: geopagg.df_alpha,
                    },
                    GeoPAGGWeightConfigEnum::RankOrder => WeightConfig::RankOrder,
                    GeoPAGGWeightConfigEnum::Balanced => WeightConfig::Balanced,
                }
            },
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
        diff_args.min_base_mean,
        misc.seed,
        &prefix,
    );
    let frame = load_dataframe(path.clone().into())?;

    let mut regex_controls = vec![];
    let mut regex_treatments = vec![];
    for label in input_args.controls {
        let regex = Regex::new(&label)?;
        regex_controls.push(regex);
    }
    for label in input_args.treatments {
        let regex = Regex::new(&label)?;
        regex_treatments.push(regex);
    }

    let mageck_results = mageck(
        &frame,
        &regex_controls,
        &regex_treatments,
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
    geopagg: GeopaggArgs,
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
            token: &misc.ntc_token,
            group_size: inc.inc_group_size,
            use_product: inc.inc_product,
            n_draws: inc.n_draws,
            fdr: misc.fdr,
        },
        GeneAggregationSelection::GeoPAGG => GeneAggregation::GeoPAGG {
            token: &misc.ntc_token,
            fdr: misc.fdr,
            weight_config: {
                match geopagg.weight_config {
                    GeoPAGGWeightConfigEnum::DropFirst => WeightConfig::DropFirst {
                        alpha: geopagg.df_alpha,
                    },
                    GeoPAGGWeightConfigEnum::RankOrder => WeightConfig::RankOrder,
                    GeoPAGGWeightConfigEnum::Balanced => WeightConfig::Balanced,
                }
            },
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
    let frame = load_dataframe(path.into())?;

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
            geopagg,
            misc,
            skip_agg,
        } => test(
            input, prefix, diff_args, agg, rra, inc, geopagg, misc, skip_agg,
        ),
        Commands::Agg {
            input,
            prefix,
            columns,
            agg,
            rra,
            inc,
            geopagg,
            misc,
        } => aggregate(input, prefix, columns, agg, rra, inc, geopagg, misc),
    }
}
