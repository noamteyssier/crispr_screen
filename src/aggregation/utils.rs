use hashbrown::HashMap;
use ndarray::Array1;

/// Converts a vector of strings to an integer representation
pub fn encode_index(genes: &Vec<String>) -> (HashMap<usize, String>, Vec<usize>)
{
    let mut total = 0usize;
    let mut map = HashMap::with_capacity(genes.len());
    let mut encoding = Vec::with_capacity(genes.len());
    for g in genes {
        match map.get(g) {
            Some(e) => {
                encoding.push(*e);
            },
            None => {
                map.insert(g, total);
                encoding.push(total);
                total += 1;
            }
        }
    }
    (map.iter().map(|(k, v)| (*v, (*k).clone())).collect(), encoding)
}

/// Select the ranks for a provided embedding. Essentially applies a filter which selects all ranks
/// for the current gene index
pub fn select_ranks(
    current_idx: usize,
    encodings: &[usize],
    ranks: &Array1<f64>) -> Array1<f64>
{
    encodings.iter()
        .zip(ranks.iter())
        .filter(|(idx, _ranks)| **idx == current_idx)
        .map(|(_, ranks)| *ranks)
        .collect()
}


#[cfg(test)]
mod testing {
    use ndarray::array;
    use super::{encode_index, select_ranks};

    #[test]
    fn test_encoding() {
        let names = vec!["g.0", "g.1", "g.0", "g.2"]
            .iter()
            .map(|x| (*x).to_string())
            .collect();
        let (encode_map, encoding) = encode_index(&names);
        assert_eq!(encoding, vec![0, 1, 0, 2]);
        assert_eq!(encode_map.get(&0), Some(&"g.0".to_string()));
        assert_eq!(encode_map.get(&1), Some(&"g.1".to_string()));
        assert_eq!(encode_map.get(&2), Some(&"g.2".to_string()));
        assert_eq!(encode_map.len(), 3);
    }

    #[test]
    fn test_selection() {
        let indices = vec![0, 1, 2, 0, 1, 2];
        let scores = array![0.5, 0., 0., 0.5, 0., 0.];
        let selection = select_ranks(0, &indices, &scores);
        assert_eq!(selection, array![0.5, 0.5]);
    }
}
