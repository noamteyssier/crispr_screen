use hashbrown::HashMap;
use ndarray::Array1;

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

pub fn alpha_rra(
    pvalues: &Array1<f64>,
    genes: &Vec<String>) 
{
    let encode = encode_index(genes);
}

#[cfg(test)]
mod testing {
    use super::encode_index;

    #[test]
    fn test_encoding() {
        let names = vec!["g.0", "g.1", "g.0", "g.2"]
            .iter()
            .map(|x| x.to_string())
            .collect();
        let encoding = encode_index(&names);
        assert_eq!(encoding, vec![0, 1, 0, 2]);
    }
}
