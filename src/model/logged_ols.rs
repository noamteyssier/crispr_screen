use std::ops::Sub;
use hashbrown::HashSet;
use ndarray::Array1;
use ndarray_rand::rand_distr::num_traits::Pow;
use super::Ols;

pub struct LoggedOls {
    kappa: f64,
    beta: f64
}
impl LoggedOls {
    pub fn fit(
        means: &Array1<f64>,
        variances: &Array1<f64>) -> Self
    {
        // subset arrays to values which won't cause numerical instability
        let (sub_means, sub_variances) = Self::subset_arrays(means, variances);

        // calculate y: log(var - mean)
        let log_variances = sub_variances.sub(&sub_means).mapv(f64::ln);

        // calculate x: log(mean)
        let log_means = sub_means.mapv(f64::ln);

        // fit ordinary least squares to log transformation
        let ols = Ols::fit(&log_means, &log_variances);
        let (kappa, beta) = (ols.alpha().exp(), ols.beta()); 

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
        variances: &Array1<f64>
        ) -> (Array1<f64>, Array1<f64>)
    {

        let idx_passing = Self::set_intersection(
            // select indices where means are greater than zero
            &Self::mask_zeros(means),

            // select indices where variances are greater than means
            &Self::mask_varied(means, variances)
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
        let mut ix = a.intersection(b).copied()
            .collect::<Vec<usize>>();
        ix.sort_unstable();
        ix
    }

    /// Return all indices where variances are greater than sample means
    fn mask_varied(
        means: &Array1<f64>,
        variances: &Array1<f64>) -> HashSet<usize>
    {
        variances
            .iter()
            .zip(means.iter())
            .enumerate()
            .filter(|(_idx, (v, m))| v > m)
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Return all indices where values are above zero
    fn mask_zeros(array: &Array1<f64>) -> HashSet<usize>
    {
        array
            .iter()
            .enumerate()
            .filter(|(_idx, x)| **x > 0.)
            .map(|(idx, _)| idx)
            .collect()
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
    use ndarray::array;
    use super::LoggedOls;

    #[test]
    fn test_mask_zeros() {
        let x = array![1., 2., 0., 3.];
        let truth = vec![0, 1, 3];
        let mask = LoggedOls::mask_zeros(&x);

        assert_eq!(mask.len(), 3);
        assert!(truth.iter().all(|x| mask.contains(x)));
    }

    #[test]
    fn test_mask_varied() {
        let x = array![1., 2., 3., 4., 5.];
        let y = array![2., 3., 0., 3., 5.];
        let truth = vec![0, 1];
        let mask = LoggedOls::mask_varied(&x, &y);

        assert_eq!(mask.len(), 2);
        assert!(truth.iter().all(|x| mask.contains(x)));
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
        let x = vec![1, 2, 4, 6].into_iter().collect();
        let y = vec![6, 4, 5, 11].into_iter().collect();
        let indices = LoggedOls::set_intersection(&x, &y);
        assert_eq!(indices, vec![4, 6]);
    }
}
