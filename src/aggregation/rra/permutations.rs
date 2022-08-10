use ndarray_rand::{rand_distr::Uniform, RandomExt};
use ndarray::Array1;
use rayon::prelude::*;
use super::{filter_alpha, robust_rank_aggregation};

fn sample_normed_ranks(
    num_ranks: usize,
    num_samples: usize) -> Array1<f64>
{
    Array1::random((num_samples,), Uniform::new(1usize, num_ranks+1))
        .iter()
        .map(|x| *x as f64 / num_ranks as f64)
        .collect()
}

pub fn run_permutations(
    num_ranks: usize,
    alpha: f64,
    npermutations: usize,
    unique_size: usize) -> Vec<f64>
{
    (0..npermutations)
        .into_par_iter()
        .map(|_| sample_normed_ranks(num_ranks, unique_size))
        .map(|choices| filter_alpha(&choices, alpha))
        .map(|filtered| robust_rank_aggregation(&filtered, unique_size))
        .collect()
}

#[cfg(test)]
mod testing {
    use ndarray_rand::rand_distr::Uniform;
    use ndarray_rand::RandomExt;
    use ndarray::Array1;
    use statrs::statistics::Statistics;
    use super::run_permutations;

    #[test]
    fn test_run(){
        let unique_size = 5;
        let num_samples = 100;
        let alpha = 0.3;
        let npermutations = 1000;
        run_permutations(num_samples, alpha, npermutations, unique_size);
    }

    #[test]
    fn test_sample_int() {
        for _ in 0..100 {
            let num_samples = 100;
            let size = 5;
            let nranks = Array1::random((size,), Uniform::new(1, num_samples))
                .iter()
                .map(|x| f64::from(*x) / f64::from(num_samples))
                .collect::<Array1<f64>>();
            assert_eq!(nranks.len(), size);
            assert!((&nranks).max() <= 1.);
            assert!(nranks.min() > 0.);
        }
    }
}
