use hashbrown::{HashSet, HashMap};
use ndarray::Array1;

/// return the indices to sort a provided array of floats
pub fn argsort<A: PartialOrd>(array: &Array1<A>) -> Array1<usize>
{
    let mut order = Array1::range(0., array.len() as f64, 1.)
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();
    order.sort_by(|a, b| 
        array[*a].partial_cmp(&array[*b]).unwrap() 
    );
    Array1::from_vec(order)
}

/// Rank an array using a temporary argsort
pub fn rank<A: PartialOrd>(array: &Array1<A>) -> Array1<usize>
{
    let order = argsort(array);
    
    argsort(&order)
}

/// return the normalized ranks of an array in place
pub fn normed_ranks(array: &Array1<f64>) -> Array1<f64>
{
    rank(array)
        .iter()
        .map(|x| (x + 1) as f64)
        .map(|x| x / array.len() as f64)
        .collect()
}

/// returns a vector of unique group sizes within the gene sets
pub fn group_sizes(array: &Vec<usize>) -> Vec<usize> {
    
    let size_map = array
        .iter()
        .fold(HashMap::new(), |mut map, x| {
            *map.entry(*x).or_insert(0) += 1usize;
            map
        });
    size_map
        .values()
        .fold(HashSet::new(), |mut set, x| {set.insert(*x); set})
        .iter().copied()
        .collect()
}

/// Returns the subarray which are below the provided alpha
pub fn filter_alpha(
    array: &Array1<f64>,
    alpha: f64) -> Array1<f64>
{
    array
        .iter()
        .filter(|x| **x < alpha).copied()
        .collect()
}

pub fn empirical_cdf(
    obs: f64,
    null: &Array1<f64>) -> f64 
{
    let size = null.len();
    let count = null
        .iter()
        .filter(|x| **x < obs)
        .count();
    (count + 1) as f64 / (size + 1) as f64

}

#[cfg(test)]
mod testing {
    use super::{argsort, normed_ranks, group_sizes, filter_alpha, empirical_cdf, rank};
    use ndarray::{Array1, array};
    use ndarray_rand::RandomExt;
    use ndarray_rand::rand_distr::Uniform;
    use statrs::statistics::Statistics;

    #[test]
    fn test_argsort() {
        let floats = array![0.3, 0.2, 0.1, 0.4];
        let order = argsort(&floats);
        assert_eq!(order, array![2, 1, 0, 3]);
    }

    #[test]
    fn test_argsort_precision() {
        let floats = array![3e-300, 2e-300, 1e-300, 4e-300];
        let order = argsort(&floats);
        assert_eq!(order, array![2, 1, 0, 3]);
    }

    #[test]
    fn test_ranks() {
        let floats = array![0.3, 0.2, 0.1, 0.4];
        let ranks = rank(&floats);
        assert_eq!(ranks, array![2, 1, 0, 3]);
    }

    #[test]
    fn test_normed_ranks() {
        let floats = array![0.3, 0.2, 0.1, 0.4];
        let truth = (array![2., 1., 0., 3.] + 1.) / 4.;
        let nranks = normed_ranks(&floats);
        assert_eq!(nranks, truth);
    }

    #[test]
    fn test_group_sizes() {
        let groups = vec![0, 0, 1, 1, 1, 2];
        let truth = vec![2, 3, 1];
        let sizes = group_sizes(&groups);
        truth
            .iter()
            .for_each(|x| assert!(sizes.contains(x)));
        assert_eq!(sizes.len(), 3);
    }

    #[test]
    fn test_filter_alpha() {
        for _ in 0..100 {
            let num_samples = 5;
            let normed_ranks = Array1::random((num_samples,), Uniform::new(0., 1.));
            let alpha = 0.3;
            let filtered = filter_alpha(&normed_ranks, alpha);
            if !filtered.is_empty() {
                assert!(filtered.max() < alpha);
            }
        }
    }

    #[test]
    fn test_empirical_cdf(){
        let size = 100;
        let obs = 0.3;
        let null = Array1::random((size,), Uniform::new(0., 1.));
        empirical_cdf(obs, &null);
    }

}
