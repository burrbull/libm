/* origin: FreeBSD /usr/src/lib/msun/src/s_cosf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 * Optimized by Bruce D. Evans.
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
use super::{k_cosf, k_sinf, rem_pio2f};
use doubled::{Doubled, AsDoubled, FromMask};

#[inline]
pub fn cosf(x: f32) -> f32 {
    /* Small multiples of pi/2 rounded to double precision. */
    let pi_1o2 = Doubled::<f32>::from_mask(0x3fc90fdb, 0xb33bbd2e);
    let pi_2o2 = Doubled::<f32>::from_mask(0x40490fdb, 0xb3bbbd2e);
    let pi_3o2 = Doubled::<f32>::from_mask(0x4096cbe4, 0xb24cde2e);
    let pi_4o2 = Doubled::<f32>::from_mask(0x40c90fdb, 0xb43bbd2e);
    let xd = x.as_doubled();

    //let x1p120 = f32::from_bits(0x7b800000); // 0x1p120f === 2 ^ 120

    let mut ix = x.to_bits();
    let sign = (ix >> 31) != 0;
    ix &= 0x7fffffff;

    if ix <= 0x3f490fda {
        /* |x| ~<= pi/4 */
        if ix < 0x39800000 {
            /* |x| < 2**-12 */
            /* raise inexact if x != 0 */
            //force_eval!(x + x1p120);
            return 1.;
        }
        return k_cosf(xd);
    }
    if ix <= 0x407b53d1 {
        /* |x| ~<= 5*pi/4 */
        if ix > 0x4016cbe3 {
            /* |x|  ~> 3*pi/4 */
            return -k_cosf(if sign { xd + pi_2o2 } else { xd - pi_2o2 });
        } else {
            if sign {
                return k_sinf(xd + pi_1o2);
            } else {
                return k_sinf(pi_1o2 - xd);
            }
        }
    }
    if ix <= 0x40e231d5 {
        /* |x| ~<= 9*pi/4 */
        if ix > 0x40afeddf {
            /* |x| ~> 7*pi/4 */
            return k_cosf(if sign { xd + pi_4o2 } else { xd - pi_4o2 });
        } else {
            if sign {
                return k_sinf(-xd - pi_3o2);
            } else {
                return k_sinf(xd - pi_3o2);
            }
        }
    }

    /* cos(Inf or NaN) is NaN */
    if ix >= 0x7f800000 {
        return f32::NAN;
    }

    /* general argument reduction needed */
    let (n, y) = rem_pio2f(x);
    let y : Doubled<f32> = y.into();
    match n & 3 {
        0 => k_cosf(y),
        1 => k_sinf(-y),
        2 => -k_cosf(y),
        _ => k_sinf(y),
    }
}
