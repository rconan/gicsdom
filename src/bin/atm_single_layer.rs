use gicsdom::ceo;
use gicsdom::ceo::Propagation;
use std::time::{Duration, Instant};


fn main() {
    let mut gmt = ceo::Gmt::new(0, None);
    gmt.build();
    let mut wfs = ceo::ShackHartmann::new(1, 48, 16, 25.5 / 48.0);
    wfs.build(8, Some(24), None);
    let mut src = wfs.new_guide_stars();
    src.build("V", vec![0.0f32], vec![0.0f32], vec![0.0f32]);

    let mut atm = ceo::Atmosphere::new();
    atm.build(
        0.15f32,
        30f32,
        1,
        vec![0f32],
        vec![1f32],
        vec![7f32],
        vec![0f32],
    );

    let n_sample: usize = 10;
    let mut wfe_rms: Vec<f32> = Vec::with_capacity(n_sample);
    let now = Instant::now();
    for k in 0..n_sample {
        src.through(&mut gmt).xpupil();
        atm.time_propagate(0.1 * k as f64, &mut src);
        wfe_rms.push(src.wfe_rms_10e(-9)[0]);
    }
    println!("Propagation in {}ms", now.elapsed().as_millis());

    let mut atm = ceo::Atmosphere::new();
    atm.raytrace_build(
        0.15f32,
        30f32,
        1,
        vec![0f32],
        vec![1f32],
        vec![7f32],
        vec![0f32],
        26f32, 1300,
        (20f32/60f32).to_radians(),
        30f32);

    let now = Instant::now();
    for k in 0..n_sample {
        src.through(&mut gmt).xpupil();
        atm.time_propagate(0.1 * k as f64, &mut src);
        wfe_rms.push(src.wfe_rms_10e(-9)[0]);
    }
    println!("Propagation in {}ms", now.elapsed().as_millis());

    println!("WFE RMS (polar-log)   [{:?}]nm", &wfe_rms[..n_sample]);
    println!("WFE RMS (ray tracing) [{:?}]nm", &wfe_rms[n_sample..]);
}
