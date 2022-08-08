use hashbrown::{HashSet, HashMap};
use ndarray::Array1;

/// Converts a vector of strings to an integer representation
pub fn encode_index(genes: &Vec<String>) -> Vec<usize>
{
    let mut total = 0usize;
    let mut map = HashMap::with_capacity(genes.len());
    let mut encoding = Vec::with_capacity(genes.len());
    for g in genes {
        match map.get(g) {
            Some(e) => {
                encoding.push(*e)
            },
            None => {
                map.insert(g, total);
                encoding.push(total);
                total += 1;
            }
        }
    }
    encoding
}

/// return the indices to sort a provided array of floats
pub fn argsort(array: &Array1<f64>) -> Vec<usize>
{
    let mut order = Array1::range(0., array.len() as f64, 1.)
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();
    order.sort_by(|a, b| 
        array[*a].partial_cmp(&array[*b]).unwrap() 
    );
    order
}

/// return the normalized ranks of an array in place
pub fn normed_ranks(array: &Array1<f64>) -> Array1<f64>
{
    argsort(array)
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
        .iter()
        .map(|x| *x)
        .collect()
}

#[cfg(test)]
mod testing {
    use super::{encode_index, argsort, normed_ranks, group_sizes};
    use ndarray::array;

    #[test]
    fn test_encoding() {
        let names = vec!["g.0", "g.1", "g.0", "g.2"]
            .iter()
            .map(|x| x.to_string())
            .collect();
        let encoding = encode_index(&names);
        assert_eq!(encoding, vec![0, 1, 0, 2]);
    }

    #[test]
    fn test_argsort() {
        let floats = array![0.3, 0.2, 0.1, 0.4];
        let order = argsort(&floats);
        assert_eq!(order, vec![2, 1, 0, 3]);
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
}
