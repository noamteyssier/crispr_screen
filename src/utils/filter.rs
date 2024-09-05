use super::logging::Logger;
use bon::builder;
use ndarray::Array2;

type FilterTuple = (Array2<f64>, Vec<String>, Vec<String>);
#[builder]
pub fn filter_low_counts(
    norm_matrix: &Array2<f64>,
    sgrna_names: &[String],
    gene_names: &[String],
    min_base: f64,
    n_controls: Option<usize>,
    logger: &Logger,
) -> FilterTuple {
    logger.filtering(min_base);

    // Calculate the mean of sgRNAs across all controls
    let sgrna_means = if let Some(n_controls) = n_controls {
        norm_matrix
            .select(ndarray::Axis(1), &(0..n_controls).collect::<Vec<usize>>())
            .mean_axis(ndarray::Axis(1))
            .expect("Failed to calculate mean of normalized matrix")
    // Calculate the mean of sgRNAs across all samples
    } else {
        norm_matrix
            .mean_axis(ndarray::Axis(1))
            .expect("Failed to calculate mean of normalized matrix")
    };
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

#[cfg(test)]
mod testing {

    use super::*;
    use ndarray::array;

    #[test]
    fn test_filtering() {
        let logger = Logger::new();
        let norm_matrix = array![
            [11.0, 12.0],
            [13.0, 14.0],
            [5.0, 6.0],
            [17.0, 18.0],
            [19.0, 20.0],
        ];
        let sgrna_names = vec![
            "sgrna1".to_string(),
            "sgrna2".to_string(),
            "sgrna3".to_string(),
            "sgrna4".to_string(),
            "sgrna5".to_string(),
        ];
        let gene_names = vec![
            "gene1".to_string(),
            "gene2".to_string(),
            "gene1".to_string(),
            "gene2".to_string(),
            "gene1".to_string(),
        ];
        let min_base = 10.0;
        let (filt_matrix, filt_sgrna_names, filt_gene_names) = filter_low_counts()
            .norm_matrix(&norm_matrix)
            .sgrna_names(&sgrna_names)
            .gene_names(&gene_names)
            .min_base(min_base)
            .logger(&logger)
            .call();

        let expected_matrix = array![[11.0, 12.0], [13.0, 14.0], [17.0, 18.0], [19.0, 20.0]];
        let expected_sgrna_names = vec!["sgrna1", "sgrna2", "sgrna4", "sgrna5"];
        let expected_gene_names = vec!["gene1", "gene2", "gene2", "gene1"];

        assert_eq!(filt_matrix, expected_matrix);
        assert_eq!(filt_sgrna_names, expected_sgrna_names);
        assert_eq!(filt_gene_names, expected_gene_names);
    }

    #[test]
    fn test_filtering_controls() {
        let logger = Logger::new();
        let norm_matrix = array![
            [11.0, 12.0],
            [13.0, 14.0],
            [5.0, 6.0],
            [17.0, 18.0],
            [19.0, 20.0],
        ];
        let sgrna_names = vec![
            "sgrna1".to_string(),
            "sgrna2".to_string(),
            "sgrna3".to_string(),
            "sgrna4".to_string(),
            "sgrna5".to_string(),
        ];
        let gene_names = vec![
            "gene1".to_string(),
            "gene2".to_string(),
            "gene1".to_string(),
            "gene2".to_string(),
            "gene1".to_string(),
        ];
        let min_base = 11.5;
        let (filt_matrix, filt_sgrna_names, filt_gene_names) = filter_low_counts()
            .norm_matrix(&norm_matrix)
            .sgrna_names(&sgrna_names)
            .gene_names(&gene_names)
            .min_base(min_base)
            .logger(&logger)
            .n_controls(1)
            .call();

        let expected_matrix = array![[13.0, 14.0], [17.0, 18.0], [19.0, 20.0]];
        let expected_sgrna_names = vec!["sgrna2", "sgrna4", "sgrna5"];
        let expected_gene_names = vec!["gene2", "gene2", "gene1"];

        assert_eq!(filt_matrix, expected_matrix);
        assert_eq!(filt_sgrna_names, expected_sgrna_names);
        assert_eq!(filt_gene_names, expected_gene_names);
    }
}
