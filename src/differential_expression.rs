use anyhow::Result;
use ndarray_rand::rand_distr::num_traits::Pow;
use statrs::distribution::{NegativeBinomial, DiscreteCDF};
use std::{fs::File, ops::{Div, Sub}};
use ndarray::{s, Array2, Axis, Array1, Zip};
use polars::prelude::{DataFrame, Series, NamedFrom, df, CsvWriter, SerWriter};
use crate::{
    utils::{parse_to_string_vec, parse_to_ndarray}, 
    math::{median_ratio_normalization, total_normalization, LoggedOLS}
};

pub enum Normalization {
    MedianRatio,
    Total
}

/// Normalize read counts using the provided method
fn normalize_counts(
    count_matrix: &Array2<f64>,
    normalization: Normalization) -> Array2<f64>
{
    match normalization{
        Normalization::MedianRatio => median_ratio_normalization(count_matrix),
        Normalization::Total => total_normalization(count_matrix)
    }
}

/// Model Mean Variance using Ordinary Least Squares Regression
fn model_mean_variance(
    normed_matrix: &Array2<f64>,
    n_controls: usize) -> Array1<f64>
{
    let model_matrix = if n_controls == 1 {
        normed_matrix.view()
    } else {
        normed_matrix.slice(s![.., ..n_controls])
    };

    let model_mean = model_matrix.mean_axis(Axis(1))
        .expect("Unexpected Empty Model Matrix");
    let model_var = model_matrix.var_axis(Axis(1), 1.);
    let logged_ols = LoggedOLS::fit(&model_mean, &model_var);
    let adj_var = logged_ols.predict(&model_mean);
    adj_var
}

fn enrichment_test(
    t_mean: f64,
    r: f64,
    p: f64,
    is_positive: bool) -> f64
{
    let distribution = NegativeBinomial::new(r, p)
        .expect("Unexpected parameterization of negative binomial");
    let pvalue = distribution.cdf(t_mean.round() as u64);
    if is_positive{
        1. - pvalue
    } else {
        pvalue
    }
}

fn enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize) -> Array1<f64>
{
    let control_means = normed_matrix
        .slice(s![.., ..n_controls])
        .mean_axis(Axis(1))
        .expect("Unexpected Empty Control Matrix")
        .map(|x| if *x == 0. {1.} else {*x});
    let treatment_means = normed_matrix
        .slice(s![.., n_controls..])
        .mean_axis(Axis(1))
        .expect("Unexpected Empty Treatment Matrix");
    
    let param_r = (&control_means.mapv(|x| x.pow(2))).div(adj_var.sub(&control_means)).mapv(|x| if x > 0. { x } else { 1. });
    let param_p = (&control_means).div(adj_var);
    let positive_lfc: Array1<bool> = (&control_means).iter().zip(treatment_means.iter()).map(|(c, t)| c < t).collect();

    Zip::from(&treatment_means)
        .and(&param_r)
        .and(&param_p)
        .and(&positive_lfc)
        .map_collect(|t_mean, r, p, is_positive| enrichment_test(*t_mean, *r, *p, *is_positive))
}

pub fn mageck(
    frame: &DataFrame,
    labels_controls: &[String],
    labels_treatments: &[String],
    normalization: Normalization
    ) -> Result<()>
{
    let columns = frame.get_column_names();
    let labels = [labels_controls, labels_treatments].concat();
    let count_matrix = parse_to_ndarray(frame, &labels)?;
    let _sgrna_names = parse_to_string_vec(frame, columns[0])?;
    let _gene_names = parse_to_string_vec(frame, columns[1])?;

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, normalization);

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&normed_matrix, labels_controls.len());

    // sgRNA Ranking (Enrichment)
    let pvalues = enrichment_testing(&normed_matrix, &adj_var, labels_controls.len());

    let mut frame = df!(
        "sgrna" => _sgrna_names,
        "gene" => _gene_names,
        "control" => normed_matrix.slice(s![.., 0]).to_vec(),
        "treatment" => normed_matrix.slice(s![.., 1]).to_vec(),
        "adj_var" => adj_var.to_vec(),
        "pvalues" => pvalues.to_vec()
    )?;

    let file = File::create("adj_var.tab")?;
    CsvWriter::new(file)
        .has_header(true)
        .with_delimiter(b'\t')
        .finish(&mut frame)?;
    
    // Gene Ranking (Aggregation)
    Ok(())
}
