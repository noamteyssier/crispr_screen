use statrs::distribution::{Beta, ContinuousCDF};
use ndarray::Array1;

/// performs robust rank aggregation
pub fn robust_rank_aggregation(
    normed_ranks: &Array1<f64>,
    num_values: usize) -> f64
{
    if normed_ranks.is_empty() { return 1. }
    let mut tmp = normed_ranks.to_vec();
    tmp.sort_floats();

    //normed_ranks
    tmp
        .iter()
        .enumerate()
        .map(|(k, n)| (k+1, num_values - k, n))
        .map(|(k, b, n)| (k as f64, b as f64, n))
        .map(|(k, b, n)| Beta::new(k, b).unwrap().cdf(*n))
        .fold(1., f64::min)
}

#[cfg(test)]
mod testing {
    use super::robust_rank_aggregation;
    use ndarray_rand::rand::{thread_rng, Rng};
    use ndarray::Array1;

    #[test]
    fn test_robust_rank(){
        let mut rng = thread_rng();
        let num = 5;
        let arr = (0..num).map(|_| rng.gen()).collect::<Array1<f64>>();
        let rra = robust_rank_aggregation(&arr, num);
        assert!(rra > 0.);
        assert!(rra <= 1.);
    }
}
