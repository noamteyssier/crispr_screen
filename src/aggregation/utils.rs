use hashbrown::HashMap;

/// Converts a vector of strings to an integer representation
pub fn encode_index(genes: &Vec<String>) -> (HashMap<usize, String>, Vec<usize>)
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
    (map.iter().map(|(k, v)| (*v, (*k).clone())).collect(), encoding)
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
        let (encode_map, encoding) = encode_index(&names);
        assert_eq!(encoding, vec![0, 1, 0, 2]);
        assert_eq!(encode_map.get(&0), Some(&"g.0".to_string()));
        assert_eq!(encode_map.get(&1), Some(&"g.1".to_string()));
        assert_eq!(encode_map.get(&2), Some(&"g.2".to_string()));
        assert_eq!(encode_map.len(), 3);
    }
}
