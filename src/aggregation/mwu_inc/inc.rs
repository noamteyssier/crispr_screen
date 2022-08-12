use hashbrown::HashMap;
use ndarray::Array1;
use crate::aggregation::utils::{encode_index, select_ranks};
use super::mann_whitney_u;

/// Validates the provided token is found one and only once in the gene set
fn validate_token(
    encode_map: &HashMap<usize, String>, 
    token: &str) -> usize
{
    let ntc_index = encode_map
        .iter()
        .filter(|(_idx, gene)| gene.contains(token))
        .map(|(idx, _gene)| *idx)
        .collect::<Vec<usize>>();
    assert_eq!(ntc_index.len(), 1, "Multiple potential genes found with provided non targeting control token");
    ntc_index[0]
}


/// Performs the Mann-Whitney U test on the dataset. This process has been called *-INC in the past
/// and so I've kept the same naming for the time being.
pub fn inc(
    pvalues: &Array1<f64>,
    genes: &Vec<String>,
    token: &str) -> (Vec<String>, Array1<f64>)
{
    let (encode_map, encode) = encode_index(genes);
    let ntc_index = validate_token(&encode_map, token);
    let ntc_values = select_ranks(ntc_index, &encode, pvalues);

    let pvalues = (0..*encode.iter().max().expect("Unexpected empty encoding"))
        .filter(|x| *x != ntc_index)
        .map(|x| select_ranks(x, &encode, pvalues))
        .map(|x| mann_whitney_u(&x, &ntc_values))
        .collect::<Array1<f64>>();

    let names = (0..*encode.iter().max().expect("Unexpected empty encoding"))
        .filter(|x| *x != ntc_index)
        .map(|curr| encode_map.get(&curr).expect("Unexpected missing index").clone())
        .collect();

    (names, pvalues)
}