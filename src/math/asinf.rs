/* origin: FreeBSD /usr/src/lib/msun/src/e_asinf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */
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

use super::consts::*;
use super::{fabsf, sqrt};
use core::f32;
use core::f64;

const PIO2: f64 = f64::consts::FRAC_PI_2;

/* coefficients for R(x^2) */
const P_S0: f32 = 1.666_658_669_7_e-01;
const P_S1: f32 = -4.274_342_209_1_e-02;
const P_S2: f32 = -8.656_363_003_0_e-03;
const Q_S1: f32 = -7.066_296_339_0_e-01;

const UF_0_5: u32 = 0x_3f00_0000;

#[inline]
fn r(z: f32) -> f32 {
    let p = z * (P_S0 + z * (P_S1 + z * P_S2));
    let q = 1. + z * Q_S1;
    p / q
}

/// Arcsine (f32)
///
/// Computes the inverse sine (arc sine) of the argument `x`.
/// Arguments to asin must be in the range -1 to 1.
/// Returns values in radians, in the range of -pi/2 to pi/2.
#[inline]
#[cfg_attr(all(test, assert_no_panic), no_panic::no_panic)]
pub fn asinf(mut x: f32) -> f32 {
    let x1p_120 = f64::from_bits(0x_3870_0000_0000_0000); // 0x1p-120 === 2 ^ (-120)

    let hx = x.to_bits();
    let ix = hx & UF_ABS;

    if ix >= UF_1 {
        /* |x| >= 1 */
        if ix == UF_1 {
            /* |x| == 1 */
            return ((x as f64) * PIO2 + x1p_120) as f32; /* asin(+-1) = +-pi/2 with inexact */
        }
        return 0. / (x - x); /* asin(|x|>1) is NaN */
    }

    if ix < UF_0_5 {
        /* |x| < 0.5 */
        /* if 0x1p-126 <= |x| < 0x1p-12, avoid raising underflow */
        if (ix < 0x_3980_0000) && (ix >= UF_MIN) {
            return x;
        }
        return x + x * r(x * x);
    }

    /* 1 > |x| >= 0.5 */
    let z = (1. - fabsf(x)) * 0.5;
    let s = sqrt(z as f64);
    x = (PIO2 - 2. * (s + s * (r(z) as f64))) as f32;
    if (hx & UF_SIGN) != 0 {
        -x
    } else {
        x
    }
}
