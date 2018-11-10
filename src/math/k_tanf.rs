/* origin: FreeBSD /usr/src/lib/msun/src/k_tan.c */
/*
 * ====================================================
 * Copyright 2004 Sun Microsystems, Inc.  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
use doubled::{Doubled, FromMask, RecPre};

/* |tan(x)/x - t(x)| < 2**-25.5 (~[-2e-08, 2e-08]). */
#[inline]
pub fn k_tanf(x: Doubled<f32>, odd: bool) -> f32 {
    let t0 = Doubled::<f32>::from_mask(0x3eaaaa6a, 0xb23e7366); // 0.333331395030791399758
    let t1 = Doubled::<f32>::from_mask(0x3e0897ea, 0xb16ccc12); // 0.133392002712976742718
    let t2 = Doubled::<f32>::from_mask(0x3d5aa649, 0xaf9e6940); // 0.0533812378445670393523
    let t3 = Doubled::<f32>::from_mask(0x3cc8ef9d, 0xb0773cc3); // 0.0245283181166547278873
    let t4 = Doubled::<f32>::from_mask(0x3b42ed70, 0xadc4c2ec); // 0.00297435743359967304927
    let t5 = Doubled::<f32>::from_mask(0x3c1b15ce, 0xad51c866); // 0.00946564784943673166728
    let z = x * x;
    /*
     * Split up the polynomial into small independent terms to give
     * opportunities for parallel evaluation.  The chosen splitting is
     * micro-optimized for Athlons (XP, X64).  It costs 2 multiplications
     * relative to Horner's method on sequential machines.
     *
     * We add the small terms from lowest degree up for efficiency on
     * non-sequential machines (the lowest degree terms tend to be ready
     * earlier).  Apart from this, we don't care about order of
     * operations, and don't need to to care since we have precision to
     * spare.  However, the chosen splitting is good for accuracy too,
     * and would give results as accurate as Horner's method if the
     * small terms were added from highest degree down.
     */
    let mut r = t4 + z * t5;
    let t = t2 + z * t3;
    let w = z * z;
    let s = z * x;
    let u = t0 + z * t1;
    r = (x + s * u) + (s * w) * (t + w * r);
    (if odd { -r.recpre() } else { r }).into()
}
