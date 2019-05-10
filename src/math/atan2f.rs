/* origin: FreeBSD /usr/src/lib/msun/src/e_atan2f.c */
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

use super::atanf;
use super::fabsf;
use crate::math::consts::*;

const PI: f32 = 3.141_592_741_0; /* 0x_4049_0fdb */
const PI_LO: f32 = -8.742_277_657_3_e-08; /* 0x_b3bb_bd2e */

/// Arctangent of y/x (f32)
///
/// Computes the inverse tangent (arc tangent) of `y/x`.
/// Produces the correct result even for angles near pi/2 or -pi/2 (that is, when `x` is near 0).
/// Returns a value in radians, in the range of -pi to pi.
#[inline]
#[cfg_attr(all(test, assert_no_panic), no_panic::no_panic)]
pub fn atan2f(y: f32, x: f32) -> f32 {
    if x.is_nan() || y.is_nan() {
        return x + y;
    }
    let mut ix = x.to_bits();
    let mut iy = y.to_bits();

    if ix == UF_1 {
        /* x = 1. */
        return atanf(y);
    }
    let m = ((iy >> 31) & 1) | ((ix >> 30) & 2); /* 2*sign(x)+sign(y) */
    ix &= UF_ABS;
    iy &= UF_ABS;

    /* when y = 0 */
    if iy == 0 {
        return match m {
            0 | 1 => y,   /* atan(+-0,+anything)=+-0 */
            2 => PI,      /* atan(+0,-anything) = pi */
            3 | _ => -PI, /* atan(-0,-anything) =-pi */
        };
    }
    /* when x = 0 */
    if ix == 0 {
        return if m & 1 != 0 { -PI / 2. } else { PI / 2. };
    }
    /* when x is INF */
    if ix == UF_INF {
        return if iy == UF_INF {
            match m {
                0 => PI / 4.,           /* atan(+INF,+INF) */
                1 => -PI / 4.,          /* atan(-INF,+INF) */
                2 => 3. * PI / 4.,      /* atan(+INF,-INF)*/
                3 | _ => -3. * PI / 4., /* atan(-INF,-INF)*/
            }
        } else {
            match m {
                0 => 0.,      /* atan(+...,+INF) */
                1 => -0.,     /* atan(-...,+INF) */
                2 => PI,      /* atan(+...,-INF) */
                3 | _ => -PI, /* atan(-...,-INF) */
            }
        };
    }
    /* |y/x| > 0x1p26 */
    if (ix + (26 << 23) < iy) || (iy == UF_INF) {
        return if m & 1 != 0 { -PI / 2. } else { PI / 2. };
    }

    /* z = atan(|y/x|) with correct underflow */
    let z = if (m & 2 != 0) && (iy + (26 << 23) < ix) {
        /*|y/x| < 0x1p-26, x < 0 */
        0.
    } else {
        atanf(fabsf(y / x))
    };
    match m {
        0 => z,                /* atan(+,+) */
        1 => -z,               /* atan(-,+) */
        2 => PI - (z - PI_LO), /* atan(+,-) */
        _ => (z - PI_LO) - PI, /* case 3 */ /* atan(-,-) */
    }
}
