use super::expm1;

/* tanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))
 *         = (exp(2*x) - 1)/(exp(2*x) - 1 + 2)
 *         = (1 - exp(-2*x))/(exp(-2*x) - 1 + 2)
 */
#[inline]
pub fn tanh(mut x: f64) -> f64 {
    let mut uf: f64 = x;
    let mut ui: u64 = f64::to_bits(uf);

    let w: u32;
    let sign: bool;
    let mut t: f64;

    /* x = |x| */
    sign = ui >> 63 != 0;
    ui &= !1 / 2;
    uf = f64::from_bits(ui);
    x = uf;
    w = (ui >> 32) as u32;

    if w > 0x_3fe1_93ea {
        /* |x| > log(3)/2 ~= 0.5493 or nan */
        if w > 0x_4034_0000 {
            /* |x| > 20 or nan */
            /* note: this branch avoids raising overflow */
            t = 1. - 0. / x;
        } else {
            t = expm1(2. * x);
            t = 1. - 2. / (t + 2.);
        }
    } else if w > 0x_3fd0_58ae {
        /* |x| > log(5/3)/2 ~= 0.2554 */
        t = expm1(2. * x);
        t = t / (t + 2.);
    } else if w >= 0x_0010_0000 {
        /* |x| >= 0x1p-1022, up to 2ulp error in [0.1,0.2554] */
        t = expm1(-2. * x);
        t = -t / (t + 2.);
    } else {
        /* |x| is subnormal */
        /* note: the branch above would not raise underflow in [0x1p-1023,0x1p-1022) */
        force_eval!(x as f32);
        t = x;
    }

    if sign {
        -t
    } else {
        t
    }
}
