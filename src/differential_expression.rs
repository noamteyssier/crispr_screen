use anyhow::Result;
use polars::prelude::DataFrame;
use crate::{
    utils::{parse_to_string_vec, parse_to_ndarray}, 
    math::{median_ratio_normalization, total_normalization}
};

pub enum Normalization {
    MedianRatio,
    Total
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
    let _normed_matrix = match normalization {
        Normalization::MedianRatio => median_ratio_normalization(&count_matrix),
        Normalization::Total => total_normalization(&count_matrix)
    };

    // Mean-Variance Modeling

    // sgRNA Ranking (Enrichment)
    
    // Gene Ranking (Aggregation)
    Ok(())
}
