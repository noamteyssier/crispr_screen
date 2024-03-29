use ndarray::{Array2, Axis};

/// Performs normalization using the average of total reads in each sample
pub fn total_normalization(matrix: &Array2<f64>) -> Array2<f64> {
    let sample_totals = matrix.sum_axis(Axis(0));

    let average_size = sample_totals.mean().expect("Unexpected Empty Input");

    let sample_factor = average_size / sample_totals;

    matrix * sample_factor
}

#[cfg(test)]
mod testing {
    use super::total_normalization;
    use ndarray::{Array2, Axis};
    use ndarray_rand::{rand_distr::Uniform, RandomExt};
    use ndarray_stats::QuantileExt;

    #[test]
    fn test_total_normalization() {
        (0..1000).for_each(|_| {
            let matrix = Array2::random((10, 4), Uniform::new(0, 5)).mapv(f64::from);
            let norm = total_normalization(&matrix);

            // take sample sums for each matrix
            let sample_sums = matrix.sum_axis(Axis(0));
            let normed_sums = norm.sum_axis(Axis(0));

            // the largest sample must be scaled down
            let largest_sample = sample_sums.argmax().unwrap();
            assert!(normed_sums[largest_sample] <= sample_sums[largest_sample]);

            // the smallest sample must be scaled up
            let smallest_sample = sample_sums.argmin().unwrap();
            assert!(normed_sums[smallest_sample] >= sample_sums[smallest_sample]);

            // matrices must be equal shape
            assert_eq!(norm.shape(), matrix.shape());
        })
    }
}
