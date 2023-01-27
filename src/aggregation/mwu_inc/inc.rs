use super::mann_whitney_u;
use crate::{
    aggregation::utils::{encode_index, select_ranks},
    utils::logging::Logger,
};
use hashbrown::HashMap;
use ndarray::Array1;

/// Validates the provided token is found one and only once in the gene set
fn validate_token(encode_map: &HashMap<usize, String>, token: &str) -> usize {
    let ntc_index = encode_map
        .iter()
        .filter(|(_idx, gene)| gene.contains(token))
        .map(|(idx, _gene)| *idx)
        .collect::<Vec<usize>>();
    assert_eq!(
        ntc_index.len(),
        1,
        "Multiple potential genes found with provided non targeting control token"
    );
    ntc_index[0]
}

/// Performs the Mann-Whitney U test on the dataset. This process has been called *-INC in the past
/// and so I've kept the same naming for the time being.
pub fn inc(
    pvalues: &Array1<f64>,
    genes: &[String],
    token: &str,
    logger: &Logger,
) -> (Vec<String>, Array1<f64>, Array1<f64>) {
    let (encode_map, encode) = encode_index(genes);
    let ntc_index = validate_token(&encode_map, token);
    let ntc_values = select_ranks(ntc_index, &encode, pvalues);
    logger.num_ntcs(ntc_values.len());

    let (scores, pvalues): (Vec<f64>, Vec<f64>) =
        (0..=*encode.iter().max().expect("Unexpected empty encoding"))
            .filter(|x| *x != ntc_index)
            .map(|x| select_ranks(x, &encode, pvalues))
            .map(|x| mann_whitney_u(&x, &ntc_values))
            .unzip();

    let names = (0..=*encode.iter().max().expect("Unexpected empty encoding"))
        .filter(|x| *x != ntc_index)
        .map(|curr| {
            encode_map
                .get(&curr)
                .expect("Unexpected missing index")
                .clone()
        })
        .collect();

    (names, Array1::from_vec(scores), Array1::from_vec(pvalues))
}

#[cfg(test)]
mod testing {
    use super::validate_token;
    use hashbrown::HashMap;

    #[test]
    fn test_validate_token() {
        let encode_map = vec![
            (0, "A"),
            (1, "B"),
            (2, "C"),
            (3, "D"),
            (4, "E"),
            (5, "F"),
            (6, "G"),
        ]
        .into_iter()
        .map(|(idx, gene)| (idx, gene.to_string()))
        .collect::<HashMap<usize, String>>();

        assert_eq!(validate_token(&encode_map, "A"), 0);
        assert_eq!(validate_token(&encode_map, "B"), 1);
        assert_eq!(validate_token(&encode_map, "C"), 2);
        assert_eq!(validate_token(&encode_map, "D"), 3);
        assert_eq!(validate_token(&encode_map, "E"), 4);
        assert_eq!(validate_token(&encode_map, "F"), 5);
        assert_eq!(validate_token(&encode_map, "G"), 6);
    }

    #[test]
    #[should_panic]
    fn test_validate_token_fail() {
        let encode_map = vec![
            (0, "A"),
            (1, "B"),
            (2, "C"),
            (3, "D"),
            (4, "E"),
            (5, "F"),
            (6, "G"),
        ]
        .into_iter()
        .map(|(idx, gene)| (idx, gene.to_string()))
        .collect::<HashMap<usize, String>>();
        validate_token(&encode_map, "H");
    }
}
