use anyhow::Result;
use std::fs::File;
use ndarray::s;
use polars::prelude::{DataFrame, Series, NamedFrom, df, CsvWriter, SerWriter};
use crate::utils::{
    parse_to_string_vec, parse_to_ndarray, 
    model_mean_variance, enrichment_testing,
    normalize_counts, Normalization
};


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
