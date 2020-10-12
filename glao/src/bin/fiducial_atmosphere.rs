use ceo::{ceo, element, Conversion, CEO};
use cirrus;
use serde_pickle as pickle;
use std::fs::File;
use std::time::Instant;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_px = 769;

    let mut gmt = ceo!(element::GMT);
    let mut src = ceo!(element::SOURCE, set_pupil_sampling = [n_px]);

    let secz = 1f32 / 30f32.to_radians().cos();
    let turb_cn2_height = [275, 425, 1250, 4000, 8000, 13000]
        .iter()
        .map(|x| *x as f32 * secz)
        .collect::<Vec<f32>>();
    let turb_cn2_xi0 = [0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751];
    let wind_speed = [5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187];
    let wind_direction = [0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462];
    let mut atm = ceo::Atmosphere::new();

    println!("Loading the atmosphere ...");
    atm.raytrace_build(
        0.14,
        25.0,
        6,
        turb_cn2_height,
        turb_cn2_xi0.to_vec(),
        wind_speed.to_vec(),
        wind_direction.to_vec(),
        25.5,
        n_px as i32,
        20f32.from_arcmin(),
        20f32,
        Some("glao_fiducial_atmosphere.bin"),
        Some(20),
    );
    println!("DONE");

    let rate = 200f64;
    let duration = 400f64;
    let n_sample = (duration * rate) as usize;
    let now = Instant::now();
    let wfe_rms = (0..n_sample)
        .map(|k| {
            atm.secs = k as f64 / rate;
            src.through(&mut gmt)
                .xpupil()
                .through(&mut atm)
                .wfe_rms_10e(-9)[0]
        })
        .collect::<Vec<f32>>();
    let et = now.elapsed().as_secs();
    println!("Elapsed time: {}s", et);

    let key = format!("GLAO/fiducial_atmosphere_wfe_rms_et{}.pkl", et);
    //    let mut file = File::create(format!("fiducial_atmosphere_wfe_rms_et{}.pkl",et)).unwrap();
    //    pickle::to_writer(&mut file, &wfe_rms, true);

    cirrus::dump("us-west-2", "gmto.modeling", &key, &wfe_rms).await?;

    Ok(())
}
