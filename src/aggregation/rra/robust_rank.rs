use ndarray::Array1;
use statrs::distribution::{Beta, ContinuousCDF};

/// Sorts the rank array in ascending order
pub fn sort_array(array: &Array1<f64>) -> Array1<f64> {
    let mut vec = array.to_vec();
    vec.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    Array1::from_vec(vec)
}

/// performs robust rank aggregation
pub fn robust_rank_aggregation(normed_ranks: &Array1<f64>, num_values: usize) -> f64 {
    if normed_ranks.is_empty() {
        return 1.;
    }

    sort_array(normed_ranks)
        .iter()
        .enumerate()
        .map(|(k, n)| (k + 1, num_values - k, n))
        .map(|(k, b, n)| (k as f64, b as f64, n))
        .map(|(k, b, n)| Beta::new(k, b).unwrap().cdf(*n))
        .fold(1., f64::min)
}

#[cfg(test)]
mod testing {
    use super::{robust_rank_aggregation, sort_array};
    use ndarray::{array, Array1};
    use ndarray_rand::rand::{thread_rng, Rng};
    use statrs::distribution::{Beta, ContinuousCDF};

    #[test]
    fn test_robust_rank() {
        let mut rng = thread_rng();
        let num = 5;
        let arr = (0..num).map(|_| rng.gen()).collect::<Array1<f64>>();
        let rra = robust_rank_aggregation(&arr, num);
        assert!(rra > 0.);
        assert!(rra <= 1.);
    }

    #[test]
    fn test_robust_rank_empty() {
        for num in 0..100 {
            let arr = array![];
            let rra = robust_rank_aggregation(&arr, num);
            assert_eq!(rra, 1.);
        }
    }

    #[test]
    fn test_robust_rank_sorted() {
        let num = 5;
        let ranks_a = array![0.5, 0.3, 0.2];
        let ranks_b = array![0.3, 0.2, 0.5];
        let rra_a = robust_rank_aggregation(&ranks_a, num);
        let rra_b = robust_rank_aggregation(&ranks_b, num);
        assert_eq!(rra_a, rra_b);
    }

    #[test]
    fn test_robust_rank_exact() {
        let num = 5;
        let ranks = array![1e-6];
        let rra = robust_rank_aggregation(&ranks, num);
        let exact = Beta::new(1., num as f64).unwrap().cdf(ranks[0]);
        assert_eq!(rra, exact);
    }

    #[test]
    fn test_sorting() {
        let arr = array![3., 2., 1.,];
        let sorted = sort_array(&arr);
        assert_eq!(sorted, array![1., 2., 3.,]);
    }
}
