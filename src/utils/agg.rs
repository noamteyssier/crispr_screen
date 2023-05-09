use hashbrown::HashMap;
use ndarray::{Array1, Axis};
use std::hash::Hash;

use super::math::weighted_mean;

/// recovers the indices of all unique values in a vector and returns a hashmap of the unique values and their indices
/// # Arguments
/// * `vec` - the vector to be searched and hashed
/// ```
pub fn unique_indices<T: Eq + Hash + Clone>(vec: &[T]) -> HashMap<T, Vec<usize>> {
    let mut map = HashMap::new();
    for (i, x) in vec.iter().enumerate() {
        map.entry(x.clone()).or_insert(Vec::new()).push(i);
    }
    map
}

/// calculates the weighted fold change for each unique value in a vector
/// # Arguments
/// * `fold_change` - the fold change values
/// * `pvalues` - the pvalues corresponding to the fold change values to be inversely weighted
/// * `map` - a hashmap of the unique values and their indices
pub fn weighted_fold_change<T: Eq + Hash + Clone>(
    fold_change: &Array1<f64>,
    pvalues: &Array1<f64>,
    map: &HashMap<T, Vec<usize>>,
) -> HashMap<T, f64> {
    map.iter()
        .map(|(k, v)| {
            let fc = fold_change.select(Axis(0), v);
            let weights = 1.0 - pvalues.select(Axis(0), v);
            let weighted_fc = weighted_mean(&fc, &weights);
            (k.clone(), weighted_fc)
        })
        .collect()
}

pub fn aggregate_fold_changes(
    gene_names: &[String],
    fold_changes: &Array1<f64>,
) -> HashMap<String, f64> {
    let index_map = unique_indices(gene_names);
    index_map
        .iter()
        .map(|(k, v)| {
            let fc = fold_changes.select(Axis(0), v);
            (k.clone(), fc.mean().unwrap())
        })
        .collect()
    
    // weighted_fold_change(fold_changes, pvalues, &map)
}

#[cfg(test)]
mod testing {
    use hashbrown::HashMap;
    use ndarray::Array1;

    #[test]
    fn test_unique_indices() {
        let vec = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut map = HashMap::new();
        map.insert(1, vec![0, 10]);
        map.insert(2, vec![1, 11]);
        map.insert(3, vec![2, 12]);
        map.insert(4, vec![3, 13]);
        map.insert(5, vec![4, 14]);
        map.insert(6, vec![5, 15]);
        map.insert(7, vec![6, 16]);
        map.insert(8, vec![7, 17]);
        map.insert(9, vec![8, 18]);
        map.insert(10, vec![9, 19]);
        assert_eq!(map, super::unique_indices(&vec));
    }

    #[test]
    fn test_weighted_fold_change() {
        use hashbrown::HashMap;
        use ndarray::Array1;

        let fc = Array1::from(vec![1., 2., 3., 4.]);
        let pvalues = Array1::from(vec![0.1, 0.2, 0.3, 0.4]);
        let mut map = HashMap::new();
        map.insert(1, vec![0, 2]);
        map.insert(2, vec![1, 3]);

        let mut expected = HashMap::new();
        expected.insert(1, ((1.0 * 0.9) + (3.0 * 0.7)) / 1.6);
        expected.insert(2, ((2.0 * 0.8) + (4.0 * 0.6)) / 1.4);

        assert_eq!(expected, super::weighted_fold_change(&fc, &pvalues, &map));
    }

    #[test]
    fn test_aggregate_fold_changes() {
        let gene_names = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let fc = Array1::from(vec![1., 2., 3., 4., 5., 6., 7., 8.]);

        let mut expected = HashMap::new();
        expected.insert("A".to_string(), ((1.0 * 0.9) + (5.0 * 0.5)) / 1.4);
        expected.insert("B".to_string(), ((2.0 * 0.8) + (6.0 * 0.4)) / 1.2);
        expected.insert("C".to_string(), ((3.0 * 0.7) + (7.0 * 0.3)) / 1.0);
        expected.insert("D".to_string(), ((4.0 * 0.6) + (8.0 * 0.2)) / 0.8);

        let result = super::aggregate_fold_changes(&gene_names, &fc);

        for (k, v) in expected.iter() {
            let v_hat = result.get(k).unwrap();
            assert!((v - v_hat).abs() < 1e-8);
        }
    }
}
