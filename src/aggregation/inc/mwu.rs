use std::ops::Div;

use ndarray::{Array1, s, concatenate, Axis};
use statrs::{statistics::{Data, RankTieBreaker, OrderStatistics}, distribution::{Normal, ContinuousCDF}};

fn merged_ranks(
    x: &Array1<f64>,
    y: &Array1<f64>) -> (Array1<f64>, Array1<f64>)
{
    let midpoint = x.len();
    let joined = concatenate(Axis(0), &[x.view(), y.view()]).unwrap().to_vec();
    let ranks = Array1::from_vec(Data::new(joined).ranks(RankTieBreaker::Average));
    (
        ranks.slice(s![..midpoint]).to_owned(),
        ranks.slice(s![midpoint..]).to_owned()
    )
}

fn u_statistic(array: &Array1<f64>) -> f64
{
    let s = array.sum();
    let n = array.len() as f64;
    s - ((n * (n+1.))/2.)
}

fn u_mean(nx: f64, ny: f64) -> f64 {
    (nx * ny) / 2.
}

fn u_std(nx: f64, ny: f64) -> f64 {
    (nx * ny * (nx + ny + 1.)).div(12.).sqrt()
}

pub fn mann_whitney_u(
    x: &Array1<f64>,
    y: &Array1<f64>) -> f64 
{
    let (ranks_x, ranks_y) = merged_ranks(x, y);

    let nx = x.len() as f64;
    let ny = y.len() as f64;

    let u_t = u_statistic(&ranks_x);
    let u_c = u_statistic(&ranks_y);
    let big_u = u_t.min(u_c);

    let m_u = u_mean(nx, ny);
    let s_u = u_std(nx, ny);

    let z_u = (big_u - m_u) / s_u;
    let pv = Normal::new(0., 1.).unwrap().cdf(z_u);

    pv
}

#[cfg(test)]
mod testing {

    use ndarray::{array, Array1};

    use super::{merged_ranks, u_statistic, mann_whitney_u};

    #[test]
    fn test_merged_ranks() {
        let x = array![1., 2., 3., 4., 5.];
        let y = array![6., 7., 8., 9., 10.,];
        let (ranks_x, ranks_y) = merged_ranks(&x, &y);
        assert_eq!(ranks_x, x);
        assert_eq!(ranks_y, y);
    }

    #[test]
    fn test_merged_ranks_ties() {
        let x = array![1., 2., 2., 4., 5.];
        let y = array![6., 7., 7., 9., 10.,];
        let (ranks_x, ranks_y) = merged_ranks(&x, &y);
        assert_eq!(ranks_x, array![1., 2.5, 2.5, 4., 5.]);
        assert_eq!(ranks_y, array![6., 7.5, 7.5, 9., 10.]);
    }

    #[test]
    fn test_u_statistic_zero() {
        let x = array![1., 2., 3., 4., 5.];
        let u = u_statistic(&x);
        assert_eq!(u, 0.);
    }

    #[test]
    fn test_u_statistic_determined() {
        let x = array![1., 5., 6., 4.];
        let u = u_statistic(&x);
        assert_eq!(u, 6.);
    }

    #[test]
    fn test_mwu() {
        let x = Array1::range(1., 6., 1.);
        let y = Array1::range(6., 11., 1.);
        println!("{:?} {:?}", x, y);
        let pv = mann_whitney_u(&x, &y);
        assert!(pv - 1.2185e-2 < 1e-6);
    }
}
