use gicsdom::ceo;
use gicsdom::ceo::Propagation;
use std::thread;
use std::time::{Duration, Instant};

fn main() {
    let mut h = vec![];

    for i in 0..3 {
        h.push(thread::spawn(move || {
            ceo::set_gpu(i as i32);
            let mut gmt = ceo::Gmt::new();
            gmt.build(0, None);
            let mut wfs = ceo::ShackHartmann::new(1, 48, 16, 25.5 / 48.0);
            wfs.build(8, Some(24), None);
            let mut src = wfs.new_guide_stars();
            src.build("V", vec![0.0f32], vec![0.0f32], vec![0.0f32]);
            let mut atm = ceo::Atmosphere::new();
            let mut atmosphere_file = format!(
                "{}",
                "/home/ubuntu/DATA/phase_screens/gmtAtmosphereL030_1579821046"
            );
            let now = Instant::now();
            atm.load_from_json(&atmosphere_file).unwrap();
            println!("Loading in {}s", now.elapsed().as_secs());
            let n_sample: usize = 10;
            let mut wfe_rms: Vec<f32> = Vec::with_capacity(n_sample);
            let now = Instant::now();
            for k in 0..n_sample {
                src.through(&mut gmt).xpupil();
                atm.time_propagate(15. * k as f64, &mut src);
                wfe_rms.push(src.wfe_rms_10e(-9)[0]);
            }
            println!("Loading in {}s", now.elapsed().as_secs());
            println!("WFE RMS [{:?}]nm",wfe_rms);
        }));
    }

    for h_ in h {
        h_.join().unwrap();
    }
}
