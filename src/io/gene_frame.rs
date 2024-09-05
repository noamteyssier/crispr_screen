use anyhow::Result;
use polars::prelude::*;
use std::{fs::File, io::BufWriter};

use crate::{
    aggregation::{AggregationResult, GeneAggregation},
    utils::{config::Configuration, logging::Logger},
};

fn build_gene_frame(results: &AggregationResult) -> Result<DataFrame, PolarsError> {
    df!(
        "gene" => results.genes(),
        "fc" => results.gene_fc().to_vec(),
        "log2fc" => results.gene_log2_fc().to_vec(),
        "score_low" => results.score_low().to_vec(),
        "pvalue_low" => results.pvalues_low().to_vec(),
        "fdr_low" => results.fdr_low().to_vec(),
        "score_high" => results.score_high().to_vec(),
        "pvalue_high" => results.pvalues_high().to_vec(),
        "fdr_high" => results.fdr_high().to_vec(),
        "pvalue" => results.pvalue().to_vec(),
        "fdr" => results.fdr().to_vec(),
        "phenotype_score" => results.phenotype_score().to_vec(),
    )
}

pub fn write_gene_frame(results: &AggregationResult, prefix: &str) -> Result<(), PolarsError> {
    let mut df = build_gene_frame(results)?;
    df.sort_in_place(["fdr"], Default::default())?;
    let writer = File::create(format!("{}.gene_results.tsv", prefix)).map(BufWriter::new)?;
    CsvWriter::new(writer)
        .with_separator(b'\t')
        .include_header(true)
        .with_quote_style(QuoteStyle::Never)
        .with_float_scientific(Some(true))
        .finish(&mut df)
}

pub fn write_hit_list(
    results: &AggregationResult,
    config: &Configuration,
    logger: &Logger,
) -> Result<(), PolarsError> {
    let mut df = build_gene_frame(results)?;
    df = match config.aggregation() {
        GeneAggregation::Inc {
            token: _,
            fdr,
            group_size: _,
            n_draws: _,
            use_product,
        } => {
            if *use_product {
                let low_mask = df
                    .column("phenotype_score")?
                    .lt(results.threshold_low().unwrap())?;
                let high_mask = df
                    .column("phenotype_score")?
                    .gt(results.threshold_high().unwrap())?;
                let mask = low_mask | high_mask;
                df.filter(&mask)
            } else {
                let low_mask = df.column("fdr_low")?.lt(*fdr)?;
                let high_mask = df.column("fdr_high")?.lt(*fdr)?;
                let mask = low_mask | high_mask;
                df.filter(&mask)
            }
        }
        GeneAggregation::AlpaRRA {
            alpha: _,
            npermutations: _,
            adjust_alpha: _,
            fdr,
        } => {
            let mask = df.column("fdr")?.lt(*fdr)?;
            df.filter(&mask)
        }
        GeneAggregation::GeoPAGG {
            token: _,
            weight_config: _,
            fdr,
            use_product: _,
        } => {
            let mask = df.column("fdr")?.lt(*fdr)?;
            df.filter(&mask)
        }
    }?;

    let num_total = df.height();
    let num_enrichments = df
        .column("log2fc")?
        .gt(0.0)?
        .iter()
        .filter(|x| x.unwrap())
        .count();
    let num_depletions = df
        .column("log2fc")?
        .lt(0.0)?
        .iter()
        .filter(|x| x.unwrap())
        .count();
    logger.hit_list(num_total, num_enrichments, num_depletions);

    df.sort_in_place(["fdr"], Default::default())?;
    let writer = File::create(format!("{}.hit_list.tsv", config.prefix())).map(BufWriter::new)?;
    CsvWriter::new(writer)
        .with_separator(b'\t')
        .include_header(true)
        .with_quote_style(QuoteStyle::Never)
        .with_float_scientific(Some(true))
        .finish(&mut df)
}
