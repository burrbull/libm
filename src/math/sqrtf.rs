/* origin: FreeBSD /usr/src/lib/msun/src/e_sqrtf.c */
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
const TINY: f32 = 1_e-30;
use super::consts::*;

pub const IF_MIN: i32 = UF_MIN as i32;

#[inline]
#[cfg_attr(all(test, assert_no_panic), no_panic::no_panic)]
pub fn sqrtf(x: f32) -> f32 {
    // On wasm32 we know that LLVM's intrinsic will compile to an optimized
    // `f32.sqrt` native instruction, so we can leverage this for both code size
    // and speed.
    llvm_intrinsically_optimized! {
        #[cfg(target_arch = "wasm32")] {
            return if x < 0. {
                ::core::f32::NAN
            } else {
                unsafe { ::core::intrinsics::sqrtf32(x) }
            }
        }
    }

    let mut ix = x.to_bits() as i32;

    /* take care of Inf and NaN */
    if (ix as u32 & UF_INF) == UF_INF {
        return x * x + x; /* sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN */
    }

    /* take care of zero */
    if ix <= 0 {
        if (ix & !(UF_SIGN as i32)) == 0 {
            return x; /* sqrt(+-0) = +-0 */
        }
        if ix < 0 {
            return (x - x) / (x - x); /* sqrt(-ve) = sNaN */
        }
    }

    /* normalize x */
    let mut m = ix >> 23;
    if m == 0 {
        /* subnormal x */
        let mut i = 0_i32;
        while ix & IF_MIN == 0 {
            ix <<= 1;
            i += 1;
        }
        m -= i - 1;
    }
    m -= 127; /* unbias exponent */
    ix = (ix & 0x_007f_ffff) | IF_MIN;
    if m & 1 == 1 {
        /* odd m, double x to make it even */
        ix += ix;
    }
    m >>= 1; /* m = [m/2] */

    /* generate sqrt(x) bit by bit */
    ix += ix;
    let mut q = 0_i32;
    let mut s = 0_i32;
    let mut r = 0x_0100_0000_u32; /* r = moving bit from right to left */

    let mut t: i32;
    while r != 0 {
        t = s + r as i32;
        if t <= ix {
            s = t + r as i32;
            ix -= t;
            q += r as i32;
        }
        ix += ix;
        r >>= 1;
    }

    /* use floating add to find out rounding direction */
    if ix != 0 {
        let mut z = 1. - TINY; /* raise inexact flag */
        if z >= 1. {
            z = 1. + TINY;
            if z > 1. {
                q += 2;
            } else {
                q += q & 1;
            }
        }
    }

    ix = (q >> 1) + 0x_3f00_0000;
    ix += m << 23;
    f32::from_bits(ix as u32)
}
