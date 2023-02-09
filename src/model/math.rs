use ndarray::Array2;

/// Determininant of a 2x2 matrix
///
/// Calculated as a * d - b * c
pub fn determinant(m: &Array2<f64>) -> f64 {
    assert_eq!(m.shape(), &[2, 2]);
    m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]]
}

/// Inverse of a 2x2 matrix
///
/// Calculated as 1 / det * [[d, -b], [-c, a]]
pub fn inverse(m: &Array2<f64>) -> Array2<f64> {
    assert_eq!(m.shape(), &[2, 2]);
    let det = determinant(m);
    let mut inv = Array2::zeros((2, 2));
    inv[[0, 0]] = m[[1, 1]] / det;
    inv[[0, 1]] = -m[[0, 1]] / det;
    inv[[1, 0]] = -m[[1, 0]] / det;
    inv[[1, 1]] = m[[0, 0]] / det;
    inv
}

#[cfg(test)]
mod testing {

    use super::*;
    use ndarray::array;

    #[test]
    fn test_determinant() {
        let m = array![[1., 2.], [3., 4.]];
        assert_eq!(determinant(&m), -2.);
    }

    #[test]
    fn test_inverse() {
        let m = array![[1., 2.], [3., 4.]];
        let inv = inverse(&m);
        let expected = array![[-2., 1.], [1.5, -0.5]];
        assert_eq!(inv, expected);
    }
}
