/* origin: FreeBSD /usr/src/lib/msun/src/k_sinf.c */
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

use doubled::{Doubled, FromMask};

/* |sin(x)/x - s(x)| < 2**-37.5 (~[-4.89e-12, 4.824e-12]). */
#[inline]
pub fn k_sinf(x: Doubled<f32>) -> f32 {
    let s1 = Doubled::<f32>::from_mask(0xbe2aaaab, 0x31b34539); // -0.166666666416265235595
    let s2 = Doubled::<f32>::from_mask(0x3c088884, 0x2f96efbb); // 0.0083333293858894631756
    let s3 = Doubled::<f32>::from_mask(0xb95007cf, 0xabb2b9dd); // -0.000198393348360966317347
    let s4 = Doubled::<f32>::from_mask(0x36366c3c, 0x29c3b46a); // 0.0000027183114939898219064
    let z = x.square();
    let w = z.square();
    let r = s3 + z * s4;
    let s = z * x;
    ((x + s * (s1 + z * s2)) + s * w * r).into()
}
