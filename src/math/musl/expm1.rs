/* origin: FreeBSD /usr/src/lib/msun/src/s_expm1.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

use crate::math::consts::*;
use core::f64;

const O_THRESHOLD: f64 = 7.097_827_128_933_839_730_96_e+02; /* 0x_4086_2E42, 0x_FEFA_39EF */
const LN2_HI: f64 = 6.931_471_803_691_238_164_90_e-01; /* 0x_3fe6_2e42, 0x_fee0_0000 */
const LN2_LO: f64 = 1.908_214_929_270_587_700_02_e-10; /* 0x_3dea_39ef, 0x_3579_3c76 */
const INVLN2: f64 = f64::consts::LOG2_E; /* 0x_3ff7_1547, 0x_652b_82fe */
/* Scaled Q's: Qn_here = 2**n * Qn_above, for R(2*z) where z = hxs = x*x/2: */
const Q1: f64 = -3.333_333_333_333_313_164_28_e-02; /* BFA11111 111110F4 */
const Q2: f64 = 1.587_301_587_254_814_601_65_e-03; /* 3F5A01A0 19FE5585 */
const Q3: f64 = -7.936_507_578_674_879_424_73_e-05; /* BF14CE19 9EAADBB7 */
const Q4: f64 = 4.008_217_827_329_362_395_52_e-06; /* 3ED0CFCA 86E65239 */
const Q5: f64 = -2.010_992_181_836_243_713_26_e-07; /* BE8AFDB7 6E09C32D */

/// Exponential, base *e*, of x-1 (f64)
///
/// Calculates the exponential of `x` and subtract 1, that is, *e* raised
/// to the power `x` minus 1 (where *e* is the base of the natural
/// system of logarithms, approximately 2.71828).
/// The result is accurate even for small values of `x`,
/// where using `exp(x)-1` would lose many significant digits.
#[inline]
pub fn expm1(mut x: f64) -> f64 {
    let hi: f64;
    let lo: f64;
    let k: i32;
    let c: f64;
    let mut t: f64;
    let mut y: f64;

    let mut ui = x.to_bits();
    let hx = ((ui >> 32) & (UF_ABS as u64)) as u32;
    let sign = (ui >> 63) as i32;

    /* filter out huge and non-finite argument */
    if hx >= 0x_4043_687A {
        /* if |x|>=56*ln2 */
        if x.is_nan() {
            return x;
        }
        if sign != 0 {
            return -1.;
        }
        if x > O_THRESHOLD {
            x *= f64::from_bits(0x_7fe0_0000_0000_0000);
            return x;
        }
    }

    /* argument reduction */
    if hx > 0x_3fd6_2e42 {
        /* if  |x| > 0.5 ln2 */
        if hx < 0x_3FF0_A2B2 {
            /* and |x| < 1.5 ln2 */
            if sign == 0 {
                hi = x - LN2_HI;
                lo = LN2_LO;
                k = 1;
            } else {
                hi = x + LN2_HI;
                lo = -LN2_LO;
                k = -1;
            }
        } else {
            k = (INVLN2 * x + if sign != 0 { -0.5 } else { 0.5 }) as i32;
            t = k as f64;
            hi = x - t * LN2_HI; /* t*ln2_hi is exact here */
            lo = t * LN2_LO;
        }
        x = hi - lo;
        c = (hi - x) - lo;
    } else if hx < 0x_3c90_0000 {
        /* |x| < 2**-54, return x */
        if hx < 0x_0010_0000 {
            force_eval!(x);
        }
        return x;
    } else {
        c = 0.;
        k = 0;
    }

    /* x is now in primary range */
    let hfx = 0.5 * x;
    let hxs = x * hfx;
    let r1 = 1. + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
    t = 3. - r1 * hfx;
    let mut e = hxs * ((r1 - t) / (6. - x * t));
    if k == 0 {
        /* c is 0 */
        return x - (x * e - hxs);
    }
    e = x * (e - c) - c;
    e -= hxs;
    /* exp(x) ~ 2^k (x_reduced - e + 1) */
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        if x < -0.25 {
            return -2. * (e - (x + 0.5));
        }
        return 1. + 2. * (x - e);
    }
    ui = ((0x3ff + k) as u64) << 52; /* 2^k */
    let twopk = f64::from_bits(ui);
    if k < 0 || k > 56 {
        /* suffice to return exp(x)-1 */
        y = x - e + 1.;
        if k == 1024 {
            y *= 2. * f64::from_bits(0x_7fe0_0000_0000_0000);
        } else {
            y *= twopk;
        }
        return y - 1.;
    }
    ui = ((0x3ff - k) as u64) << 52; /* 2^-k */
    let uf = f64::from_bits(ui);
    if k < 20 {
        y = (x - e + (1. - uf)) * twopk;
    } else {
        y = (x - (e + uf) + 1.) * twopk;
    }
    y
}

#[cfg(test)]
mod tests {
    #[test]
    fn sanity_check() {
        assert_eq!(super::expm1(1.1), 2.004_166_023_946_433_4);
    }
}
