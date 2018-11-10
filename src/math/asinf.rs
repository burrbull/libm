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

use core::f32;

use super::fabsf::fabsf;
use doubled::{Doubled, FromMask, Scale};
use math::sqrtf_as_doubled;

const PIO2: f32 = 1.570796326794896558e+00;

/* coefficients for R(x^2) */
const P_S0: f32 = 1.6666586697e-01;
const P_S1: f32 = -4.2743422091e-02;
const P_S2: f32 = -8.6563630030e-03;
const Q_S1: f32 = -7.0662963390e-01;

#[inline]
fn r(z: f32) -> f32 {
    let p = z * (P_S0 + z * (P_S1 + z * P_S2));
    let q = 1. + z * Q_S1;
    p / q
}

#[inline]
pub fn asinf(mut x: f32) -> f32 {
    let pio2 = Doubled::<f32>::from_mask(0x3fc90fdb, 0xb33bbd2e);

    let hx = x.to_bits();
    let ix = hx & 0x7fffffff;
    let sign = hx & 0x80000000;

    if ix >= 0x3f800000 {
        /* |x| >= 1 */
        if ix == 0x3f800000 {
            /* |x| == 1 */
            f32::from_bits(PIO2.to_bits() | sign) /* asin(+-1) = +-pi/2 with inexact */
        } else {
            f32::NAN /* asin(|x|>1) is NaN */
        }
    } else if ix < 0x3f000000 {
        /* |x| < 0.5 */
        /* if 0x1p-126 <= |x| < 0x1p-12, avoid raising underflow */
        if (ix < 0x39800000) && (ix >= 0x00800000) {
            x
        } else {
            x + x * r(x * x)
        }
    } else {
        /* 1 > |x| >= 0.5 */
        let z = (1. - fabsf(x)) * 0.5;
        let s = sqrtf_as_doubled(z);
        x = (pio2 - (s + s * Doubled::from(r(z))).scale(2.)).into();
        if (hx >> 31) != 0 {
            -x
        } else {
            x
        }
    }
}
