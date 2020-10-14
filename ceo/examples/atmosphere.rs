use ceo::{ceo, element::*, CEO};

fn main() {
    let mut gmt = ceo!(GMT);
    let mut src = ceo!(SOURCE);
    let mut atm = ceo!(ATMOSPHERE);
    let dt = 10_f64;
    for k in 0..10 {
        atm.secs = k as f64 * dt;
        src.through(&mut gmt).xpupil().through(&mut atm);
        println!(
            "T: {:02}s -> WFE RMS: {:.0}nm",
            atm.secs,
            src.wfe_rms_10e(-9)[0]
        );
    }
}
