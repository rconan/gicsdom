use ceo::{ceo, element, Conversion, CEO};
use cirrus;
use glao::glao_sys::ScienceField;
use serde::{Deserialize, Serialize};
use std::time::Instant;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Results {
    time: Vec<f64>,
    wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>,
    atm_fwhm_x: f64,
    fwhm: Vec<f64>,
    pssn: Vec<f32>,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_px = 769;

    let mut gmt = ceo!(element::GMT);
    //    let mut src = ceo!(element::SOURCE, set_pupil_sampling = [n_px]);
    let mut science = ScienceField::delaunay_21("Vs", n_px, None);
    //let mut science = ScienceField::on_axis("Vs", n_px, None);
    science.build();

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
        science.pssn.r0(),
        science.pssn.oscale,
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
    let mut t: Vec<f64> = vec![];
    let wfe = (0..n_sample)
        .map(|k| {
            atm.secs = k as f64 / rate;
            t.push(atm.secs);
            let wfe_rms = science
                .src
                .through(&mut gmt)
                .xpupil()
                .through(&mut atm)
                .through(&mut science.pssn)
                .wfe_rms_10e(-9);
            (
                wfe_rms,
                science.src.segment_wfe_rms_10e(-9),
                science.src.segment_piston_10e(-9),
            )
        })
        .collect::<Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>>();
    let et = now.elapsed().as_secs();
    println!("Elapsed time: {}s", et);

    let key = "GLAO/fiducial_atmosphere_wfe.pkl";
    //    let mut file = File::create(format!("fiducial_atmosphere_wfe_rms_et{}.pkl",et)).unwrap();
    //    pickle::to_writer(&mut file, &wfe_rms, true);

    //    cirrus::dump("us-west-2", "gmto.modeling", &key, &wfe_rms).await?;
    let (atm_fwhm_x, fwhm) = science.glao_wrap_up();
    let results = Results {
        time: t,
        wfe,
        atm_fwhm_x,
        fwhm,
        pssn: science.pssn.peek().estimates.clone(),
    };
    println!("Uploading results to: s3://gmto.modeling/{}", key);
    cirrus::dump("us-west-2", "gmto.modeling", &key, &results).await?;

    Ok(())
}
