use crate::math::consts::*;

pub fn modff(x: f32) -> (f32, f32) {
    let rv2: f32;
    let mut u: u32 = x.to_bits();
    let mask: u32;
    let e = ((u >> 23 & 0xff) as isize) - 0x7f;

    /* no fractional part */
    if e >= 23 {
        rv2 = x;
        if e == 0x80 && (u << 9) != 0 {
            /* nan */
            return (x, rv2);
        }
        u &= UF_SIGN;
        return (f32::from_bits(u), rv2);
    }
    /* no integral part */
    if e < 0 {
        u &= UF_SIGN;
        rv2 = f32::from_bits(u);
        return (x, rv2);
    }

    mask = 0x_007f_ffff >> e;
    if (u & mask) == 0 {
        rv2 = x;
        u &= UF_SIGN;
        return (f32::from_bits(u), rv2);
    }
    u &= !mask;
    rv2 = f32::from_bits(u);
    (x - rv2, rv2)
}
