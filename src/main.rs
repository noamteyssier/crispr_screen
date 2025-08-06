use adjustp::Procedure;
use anyhow::Result;
use bon::builder;
use clap::Parser;
use cli::{
    Cli, Commands, DiffAbundanceArgs, GeopaggArgs, IncArgs, InputArgs, MiscArgs, RraArgs,
    SgrnaColumns,
};
use geopagg::WeightConfig;
use log::LevelFilter;
use run_aggregation::run_aggregation;
use std::path::Path;

pub mod aggregation;
pub mod cli;
pub mod differential_expression;
pub mod enrich;
pub mod io;
pub mod model;
pub mod norm;
pub mod resample;
pub mod run_aggregation;
pub mod utils;

use aggregation::{GeneAggregation, GeneAggregationSelection, GeoPAGGWeightConfigEnum};
use differential_expression::mageck;
use io::{build_regex_set, load_dataframe};
use resample::resample;
use utils::{config::Configuration, logging::Logger, Adjustment};

#[builder]
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
        GeneAggregationSelection::GeoPAGG => {
            let token = if misc.ntc_token.is_empty() | geopagg.use_all {
                None
            } else {
                Some(misc.ntc_token.as_str())
            };
            GeneAggregation::GeoPAGG {
                token,
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
                use_product: geopagg.use_product,
                zscore_threshold: geopagg.zscore_threshold,
            }
        }
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

    let config = Configuration::builder()
        .normalization(diff_args.norm)
        .aggregation(agg)
        .correction(correction)
        .model_choice(diff_args.model_choice)
        .min_base_mean(diff_args.min_base_mean)
        .strategy(diff_args.strategy)
        .seed(misc.seed)
        .prefix(&prefix)
        .build();
    let frame = load_dataframe(path.clone().into())?;

    let regex_controls = build_regex_set(&input_args.controls)?;
    let regex_treatments = build_regex_set(&input_args.treatments)?;

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

#[builder]
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
        GeneAggregationSelection::GeoPAGG => {
            let token = if misc.ntc_token.is_empty() | geopagg.use_all {
                None
            } else {
                Some(misc.ntc_token.as_str())
            };
            GeneAggregation::GeoPAGG {
                token,
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
                use_product: geopagg.use_product,
                zscore_threshold: geopagg.zscore_threshold,
            }
        }
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

    let config = Configuration::builder()
        .aggregation(agg)
        .correction(correction)
        .seed(misc.seed)
        .prefix(&prefix)
        .build();
    let frame = load_dataframe(path.into())?;

    run_aggregation(&frame, columns, &config, &logger)
}

fn main() -> Result<()> {
    let args = Cli::parse();

    env_logger::builder()
        .format_timestamp_secs()
        .filter_level(LevelFilter::Info)
        .parse_env("CRISPR_SCREEN_LOG")
        .init();

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
        } => test()
            .input_args(input)
            .prefix(prefix)
            .diff_args(diff_args)
            .agg(agg)
            .rra(rra)
            .inc(inc)
            .geopagg(geopagg)
            .misc(misc)
            .skip_agg(skip_agg)
            .call(),
        Commands::Agg {
            input,
            prefix,
            columns,
            agg,
            rra,
            inc,
            geopagg,
            misc,
        } => aggregate()
            .input(input)
            .prefix(prefix)
            .columns(columns)
            .agg(agg)
            .rra(rra)
            .inc(inc)
            .geopagg(geopagg)
            .misc(misc)
            .call(),
        Commands::Resample {
            input,
            output,
            n_resamples,
            depth,
            depth_samples,
            seed,
            samples,
            quiet,
        } => resample()
            .input(input)
            .maybe_path(output)
            .n_resamples(n_resamples)
            .seed(seed)
            .maybe_depth(depth)
            .maybe_depth_samples(depth_samples)
            .samples(samples)
            .quiet(quiet)
            .call(),
    }
}
