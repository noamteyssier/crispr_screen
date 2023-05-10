use crate::utils::logging::Logger;
use hashbrown::HashSet;
use ndarray::Array1;

/// Return all indices where values are above zero
pub fn mask_zeros(array: &Array1<f64>, logger: &Logger) -> HashSet<usize> {
    let mask = array
        .iter()
        .enumerate()
        .filter(|(_idx, x)| **x > 0.)
        .map(|(idx, _)| idx)
        .collect::<HashSet<usize>>();
    logger.num_zeros(array.len() - mask.len());
    mask
}

/// Select from vector where indices are in the mask
pub fn select_from_mask<T: Clone>(array: &[T], mask: &[usize]) -> Vec<T> {
    mask.iter().map(|x| array[*x].clone()).collect::<Vec<T>>()
}

/// Select from array where indices are in the mask
pub fn select_from_mask_array<T: Clone>(array: &Array1<T>, mask: &[usize]) -> Array1<T> {
    mask.iter()
        .map(|x| array[*x].clone())
        .collect::<Array1<T>>()
}

/// Filter `sgRNAs` with zero counts in both samples
pub fn filter_zeros(
    base_means: &Array1<f64>,
    gene_names: &[String],
    sgrna_pvalues_low: &Array1<f64>,
    sgrna_pvalues_high: &Array1<f64>,
    sgrna_log_fold_change: &Array1<f64>,
    logger: &Logger,
) -> (Vec<String>, Array1<f64>, Array1<f64>, Array1<f64>) {
    let passing_indices = mask_zeros(base_means, logger);
    let mut sorted_indices = passing_indices.iter().copied().collect::<Vec<usize>>();
    sorted_indices.sort_unstable();

    let passing_gene_names = select_from_mask(gene_names, &sorted_indices);
    let passing_sgrna_pvalues_low = select_from_mask_array(sgrna_pvalues_low, &sorted_indices);
    let passing_sgrna_pvalues_high = select_from_mask_array(sgrna_pvalues_high, &sorted_indices);
    let passing_log_fold_change = select_from_mask_array(sgrna_log_fold_change, &sorted_indices);

    (
        passing_gene_names,
        passing_sgrna_pvalues_low,
        passing_sgrna_pvalues_high,
        passing_log_fold_change,
    )
}

/// Calculates an empirical alpha threshold for RRA
fn calculate_empirical_alpha(pvalue_arr: &Array1<f64>, alpha: f64) -> f64 {
    pvalue_arr
        .map(|x| if *x < alpha { 1.0 } else { 0.0 })
        .mean()
        .expect("Error calculating mean in empirical alpha")
}

/// Sets the alpha threshold empirically or returns the untouched alpha
pub fn set_alpha_threshold(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    alpha: f64,
    adjust_alpha: bool,
) -> (f64, f64) {
    if adjust_alpha {
        (
            calculate_empirical_alpha(pvalue_low, alpha),
            calculate_empirical_alpha(pvalue_high, alpha),
        )
    } else {
        (alpha, alpha)
    }
}

/// Calculates the number of unique values in a vector
pub fn num_unique(names: &[String]) -> usize {
    names.iter().collect::<HashSet<_>>().len()
}

#[cfg(test)]
mod testing {
    use super::{calculate_empirical_alpha, filter_zeros, mask_zeros};
    use crate::{
        aggregation::utils::{select_from_mask, select_from_mask_array, set_alpha_threshold},
        utils::logging::Logger,
    };
    use ndarray::{array, Array1, Array2, Axis};
    use ndarray_rand::{
        rand_distr::{Binomial, Uniform},
        RandomExt,
    };

    #[test]
    fn test_empirical_alpha() {
        let array = Array1::from_vec(vec![0.1, 0.2, 0.25, 0.5, 0.5]);
        let alpha = 0.3;
        let emp = calculate_empirical_alpha(&array, alpha);
        assert_eq!(emp, 3.0 / 5.0);
    }

    #[test]
    fn test_select_from_mask() {
        let array = vec![0, 1, 2, 3, 4, 5];
        let mask = vec![0, 2, 4];
        let selected = select_from_mask(&array, &mask);
        assert_eq!(selected, vec![0, 2, 4]);
    }

    #[test]
    fn test_select_from_mask_array() {
        let array = Array1::from_vec(vec![0, 1, 2, 3, 4, 5]);
        let mask = vec![0, 2, 4];
        let selected = select_from_mask_array(&array, &mask);
        assert_eq!(selected, array![0, 2, 4]);
    }

    #[test]
    fn test_mask_zeros() {
        let array = Array1::from_vec(vec![0., 1., 0., 1., 0.]);
        let logger = Logger::new();
        let mask = mask_zeros(&array, &logger);
        assert_eq!(mask.len(), 2);
        assert!(mask.contains(&1));
        assert!(mask.contains(&3));
    }

    #[test]
    fn test_filter_zeros() {
        let logger = Logger::new();
        let array = Array2::random((100, 2), Binomial::new(1, 0.2).unwrap()).mapv(|x| x as f64);
        let means = array.mean_axis(Axis(1)).unwrap();
        let nonzero = mask_zeros(&means, &logger);
        let gene_names = (0..100)
            .map(|x| format!("gene_{x}"))
            .collect::<Vec<String>>();
        let p_low = Array1::random(100, Uniform::new(0.0, 1.0));
        let p_high = Array1::random(100, Uniform::new(0.0, 1.0));
        let logfc = Array1::random(100, Uniform::new(0.0, 1.0));
        let (pgn, ppl, pph, ppf) =
            filter_zeros(&means, &gene_names, &p_low, &p_high, &logfc, &logger);

        assert_eq!(pgn.len(), nonzero.len());
        assert_eq!(ppl.len(), nonzero.len());
        assert_eq!(pph.len(), nonzero.len());
        assert_eq!(ppf.len(), nonzero.len());
    }

    #[test]
    fn test_filter_zeros_empty() {
        let logger = Logger::new();
        let array = Array2::random((100, 2), Binomial::new(1, 0.0).unwrap()).mapv(|x| x as f64);
        let means = array.mean_axis(Axis(1)).unwrap();
        let nonzero = mask_zeros(&means, &logger);
        let gene_names = (0..100)
            .map(|x| format!("gene_{x}"))
            .collect::<Vec<String>>();
        let p_low = Array1::random(100, Uniform::new(0.0, 1.0));
        let p_high = Array1::random(100, Uniform::new(0.0, 1.0));
        let logfc = Array1::random(100, Uniform::new(0.0, 1.0));
        let (pgn, ppl, pph, ppf) =
            filter_zeros(&means, &gene_names, &p_low, &p_high, &logfc, &logger);

        assert_eq!(pgn.len(), nonzero.len());
        assert_eq!(ppl.len(), nonzero.len());
        assert_eq!(pph.len(), nonzero.len());
        assert_eq!(ppf.len(), nonzero.len());
    }

    #[test]
    fn test_filter_zeros_all() {
        let logger = Logger::new();
        let array = Array2::random((100, 2), Binomial::new(1, 1.0).unwrap()).mapv(|x| x as f64);
        let means = array.mean_axis(Axis(1)).unwrap();
        let nonzero = mask_zeros(&means, &logger);
        let gene_names = (0..100)
            .map(|x| format!("gene_{x}"))
            .collect::<Vec<String>>();
        let p_low = Array1::random(100, Uniform::new(0.0, 1.0));
        let p_high = Array1::random(100, Uniform::new(0.0, 1.0));
        let logfc = Array1::random(100, Uniform::new(0.0, 1.0));
        let (pgn, ppl, pph, ppf) =
            filter_zeros(&means, &gene_names, &p_low, &p_high, &logfc, &logger);

        assert_eq!(pgn.len(), nonzero.len());
        assert_eq!(ppl.len(), nonzero.len());
        assert_eq!(pph.len(), nonzero.len());
        assert_eq!(ppf.len(), nonzero.len());
    }

    #[test]
    fn test_set_alpha_threshold() {
        let alpha = 0.25;
        let empirical_alpha = 0.4;
        let pvalue_low = Array1::from_vec(vec![0.1, 0.2, 0.25, 0.5, 0.5]);
        let pvalue_high = Array1::from_vec(vec![0.1, 0.2, 0.25, 0.5, 0.5]);

        let (t_low, t_high) = set_alpha_threshold(&pvalue_low, &pvalue_high, alpha, true);
        assert_eq!(t_low, empirical_alpha);
        assert_eq!(t_high, empirical_alpha);
    }

    #[test]
    fn test_set_alpha_threshold_no_adjust() {
        let alpha = 0.25;
        let pvalue_low = Array1::from_vec(vec![0.1, 0.2, 0.25, 0.5, 0.5]);
        let pvalue_high = Array1::from_vec(vec![0.1, 0.2, 0.25, 0.5, 0.5]);

        let (t_low, t_high) = set_alpha_threshold(&pvalue_low, &pvalue_high, alpha, false);
        assert_eq!(t_low, alpha);
        assert_eq!(t_high, alpha);
    }
}
