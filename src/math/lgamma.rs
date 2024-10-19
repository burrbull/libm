use super::lgamma_r;

#[cfg_attr(all(test, assert_no_panic), no_panic::no_panic)]
pub const fn lgamma(x: f64) -> f64 {
    lgamma_r(x).0
}
