use ceo::{ceo, ATMOSPHERE, Builder};
use serde_pickle as pkl;
use std::fs::File;

fn main() {
    let mut gmt = ceo!(GMT);
    let mut src = ceo!(SOURCE);

    let mut atm_1 = ATMOSPHERE::new()
        .set_single_turbulence_layer(0f32, None, None)
        .build();
    let mut atm_2 = ATMOSPHERE::new()
        .set_single_turbulence_layer(0f32, None, None)
        .set_ray_tracing(25.5, 512, 0., 1., Some("atm_2.bin".to_owned()), None)
        .build();

    let dump = |data: &Vec<f32>, filename: &str| {
        let mut file = File::create(filename).unwrap();
        pkl::to_writer(&mut file, data, true).unwrap();
    };
    dump(&(src.through(&mut gmt).xpupil().through(&mut atm_1).phase()),"atm_1.pkl");
    dump(&(src.through(&mut gmt).xpupil().through(&mut atm_2).phase()),"atm_2.pkl");
    /*
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
    */
}
