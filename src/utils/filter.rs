use ndarray::Array2;

use super::logging::Logger;

type FilterTuple = (Array2<f64>, Vec<String>, Vec<String>);

pub fn filter_low_counts(
    norm_matrix: &Array2<f64>,
    sgrna_names: &[String],
    gene_names: &[String],
    min_base: f64,
    logger: &Logger,
) -> FilterTuple {
    logger.filtering(min_base);
    let sgrna_means = norm_matrix
        .mean_axis(ndarray::Axis(1))
        .expect("Failed to calculate mean of normalized matrix");
    let mask = sgrna_means
        .iter()
        .enumerate()
        .filter(|(_idx, v)| **v >= min_base)
        .map(|(idx, _)| idx)
        .collect::<Vec<usize>>();
    let filt_matrix = norm_matrix.select(ndarray::Axis(0), &mask);
    let filt_sgrna_names = mask.iter().map(|idx| sgrna_names[*idx].clone()).collect();
    let filt_gene_names = mask.iter().map(|idx| gene_names[*idx].clone()).collect();
    let num_filtered = sgrna_names.len() - mask.len();
    logger.num_filtered(num_filtered);
    (filt_matrix, filt_sgrna_names, filt_gene_names)
}
