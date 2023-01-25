use std::ops::Sub;
use hashbrown::HashSet;
use ndarray::Array1;
use ndarray_rand::rand_distr::num_traits::Pow;
use crate::utils::{math::zscore_transform, logging::Logger};
use super::{Ols, ModelChoice, Wols};

pub struct LoggedOls {
    kappa: f64,
    beta: f64
}
impl LoggedOls {
    pub fn fit(
        means: &Array1<f64>,
        variances: &Array1<f64>,
        model_choice: &ModelChoice,
        logger: &Logger) -> Self
    {
        logger.start_mean_variance();

        // subset arrays to values which won't cause numerical instability
        let (sub_means, sub_variances) = Self::subset_arrays(means, variances, logger);

        // calculate y: log(var - mean)
        let log_variances = sub_variances.sub(&sub_means).mapv(|x| x + 1.).mapv(f64::ln);

        // calculate x: log(mean)
        let log_means = sub_means.mapv(|x| x + 1.).mapv(f64::ln);

        let (kappa, beta) = match model_choice {
            ModelChoice::Ols => {
                let ols = Ols::fit(&log_means, &log_variances);
                (ols.alpha().exp(), ols.beta())
            },
            ModelChoice::Wols => {
                // use non-logged means as the weights 
                let wols = Wols::fit(&log_means, &log_variances, &sub_means);
                (wols.alpha().exp(), wols.beta())
            }
        };

        logger.ols_parameters(model_choice, kappa, beta);
        Self { kappa, beta }
    }

    /// Calculates the adjusted variance from the model fit parameters
    pub fn predict(&self, means: &Array1<f64>) -> Array1<f64>
    {
        // map adjusted variance formula as:
        // adj_var = means + (kappa * (means ** beta))
        let adj_var = means + (self.kappa * (means.mapv(|x| x.pow(self.beta))));

        // replace all zeros with the nonzero minimum
        Self::replace_zeros_with_min(&adj_var)
    }

    /// Subset arrays to those that will not cause numerical instability
    fn subset_arrays(
        means: &Array1<f64>,
        variances: &Array1<f64>,
        logger: &Logger) -> (Array1<f64>, Array1<f64>)
    {

        let idx_passing = Self::set_intersection(

            // select indices where means are under 4 std of global mean
            &Self::mask_outliers(means, logger),

            // select indices where variances are greater than means
            &Self::mask_varied(means, variances, logger)
            );

        (
            idx_passing.iter().map(|idx| means[*idx]).collect(),
            idx_passing.iter().map(|idx| variances[*idx]).collect()
        )
    }

    /// Return all unique indices 
    fn set_intersection(
        a: &HashSet<usize>, 
        b: &HashSet<usize>) -> Vec<usize>
    {
        let mut ix = (a & b).into_iter().collect::<Vec<usize>>();
        ix.sort_unstable();
        ix
    }

    /// Return all indices where means are less than 4 standard deviations away from the global
    /// mean
    fn mask_outliers(means: &Array1<f64>, logger: &Logger) -> HashSet<usize>
    {
        let mask = zscore_transform(means)
            .iter()
            .enumerate()
            .filter(|(_idx, x)| **x < 4.)
            .map(|(idx, _)| idx)
            .collect::<HashSet<usize>>();
        logger.num_outliers(means.len() - mask.len());
        mask
    }

    /// Return all indices where variances are greater than sample means
    fn mask_varied(
        means: &Array1<f64>,
        variances: &Array1<f64>,
        logger: &Logger) -> HashSet<usize>
    {
        let mask = variances
            .iter()
            .zip(means.iter())
            .enumerate()
            .filter(|(_idx, (v, m))| v > m)
            .map(|(idx, _)| idx)
            .collect::<HashSet<usize>>();
        logger.num_varied(means.len() - mask.len());
        mask
    }

    /// Calculates the minimum value in an array
    fn min_nonzero(array: &Array1<f64>) -> Option<&f64>
    {
        array
            .iter()
            .filter(|x| **x > 0.)
            .reduce(|a, b| if b < a { b } else { a })
    }

    /// Replace all zero elements with the minumum nonzero array
    fn replace_zeros_with_min(array: &Array1<f64>) -> Array1<f64>
    {
        if let Some(min) = Self::min_nonzero(array) {
            array.mapv(|x| if x > 0. { x } else { *min })
        }
        else {
            array.to_owned()
        }
    }
}

#[cfg(test)]
mod testing {
    use ndarray::{array, Array1, concatenate, Axis};
    use ndarray_rand::{RandomExt, rand_distr::{Normal, Uniform}};
    use crate::utils::{math::zscore_transform, logging::Logger};

    use super::LoggedOls;

    #[test]
    fn test_mask_varied() {
        let x = array![1., 2., 3., 4., 5.];
        let y = array![2., 3., 0., 3., 5.];
        let truth = vec![0, 1];
        let logger = Logger::new();
        let mask = LoggedOls::mask_varied(&x, &y, &logger);

        assert_eq!(mask.len(), 2);
        assert!(truth.iter().all(|x| mask.contains(x)));
    }

    #[test]
    fn test_mask_zscore() {
        let x = Array1::random(1000, Normal::new(0., 1.).unwrap());
        let y = Array1::random(5, Uniform::new(100., 150.));
        let merged = concatenate(Axis(0), &[x.view(), y.view()]).unwrap();
        let znorm = zscore_transform(&merged);
        let logger = Logger::new();
        let passing = LoggedOls::mask_outliers(&znorm, &logger);
        assert_eq!(passing.len(), 1000);
        (1000..1005).for_each(|x| {
            assert!(!passing.contains(&x));
        });
    }

    #[test]
    fn test_replace_nonzero() {
        let x = array![1., 2., 0., 3.];
        let y = LoggedOls::replace_zeros_with_min(&x);
        assert_eq!(y, array![1., 2., 1., 3.]);
    }

    #[test]
    fn test_min_nonzero_some() {
        let x = array![1., 2., 0., 3.];
        let m = LoggedOls::min_nonzero(&x);
        assert_eq!(m, Some(&1.));
    }


    #[test]
    fn test_min_nonzero_none() {
        let x = array![0., 0., 0.];
        let m = LoggedOls::min_nonzero(&x);
        assert_eq!(m, None);
    }


    #[test]
    fn test_sorted_intersection() {
        let x = vec![1, 2, 4, 6, 12].into_iter().collect();
        let y = vec![6, 4, 5, 12, 11].into_iter().collect();
        let indices = LoggedOls::set_intersection(&x, &y);
        assert_eq!(indices, vec![4, 6, 12]);
    }
}
