/* origin: FreeBSD /usr/src/lib/msun/src/k_cosf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 * Debugged and optimized by Bruce D. Evans.
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

/* |cos(x) - c(x)| < 2**-34.1 (~[-5.37e-11, 5.295e-11]). */
#[inline]
pub fn k_cosf(x: Doubled<f32>) -> f32 {
    let c0 = Doubled::<f32>::from_mask(0xbf000000, 0x313ce860); // -0.499999997251031003120
    let c1 = Doubled::<f32>::from_mask(0x3d2aaa9f, 0x2f029d21); // 0.0416666233237390631894
    let c2 = Doubled::<f32>::from_mask(0xbab6043f, 0xae00f1e2); // -0.00138867637746099294692
    let c3 = Doubled::<f32>::from_mask(0x37cc9a17, 0x296e5069); // 0.0000243904487962774090654
    let z = x.square();
    let w = z.square();
    let r = c2 + z * c3;
    (((1.0 + z * c0) + w * c1) + (w * z) * r).into()
}
