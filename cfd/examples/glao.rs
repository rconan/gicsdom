use ceo::{
    calibrations, ceo, element::*, pssn::AtmosphereTelescopeError as ATE,
    shackhartmann::Geometric as WFS_TYPE, CEOInit, Calibration, Conversion, CEO,
};
use cfd;
use indicatif::ProgressBar;
//use serde_pickle as pickle;
//use std::fs::File;
use std::time::Instant;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_px = 769;
    let n_kl = 70;
    let gmt_blueprint = CEO::<GMT>::new().set_m2_n_mode(n_kl);
    let wfs_blueprint = CEO::<SH48>::new();
    let gs_blueprint = wfs_blueprint.guide_stars().set_on_ring(6f32.from_arcmin());

    let mut gmt2wfs = Calibration::new(&gmt_blueprint, &gs_blueprint, &wfs_blueprint);
    let mirror = vec![calibrations::Mirror::M2MODES];
    let segments = vec![vec![calibrations::Segment::Modes(1e-6, 1..n_kl)]; 7];
    let now = Instant::now();
    let wfs_intensity_threshold = 0.8;
    gmt2wfs.calibrate(mirror, segments, Some(wfs_intensity_threshold));
    println!(
        "GTM 2 WFS calibration [{}x{}] in {}s",
        gmt2wfs.n_data,
        gmt2wfs.n_mode,
        now.elapsed().as_secs()
    );
    let poke_sum = gmt2wfs.poke.from_dev().iter().sum::<f32>();
    println!("Poke sum: {}", poke_sum);

    let mut gmt = gmt_blueprint.build();
    let mut wfs = wfs_blueprint.build::<WFS_TYPE>();
    let mut gs = gs_blueprint.build();

    gs.through(&mut gmt).xpupil();
    println!("GS WFE RMS: {}nm", gs.wfe_rms_10e(-9)[0]);
    wfs.calibrate(&mut gs, wfs_intensity_threshold);
    println!("# valid lenslet: {}", wfs.n_valid_lenslet());

    let mut src = ceo!(SOURCE, set_pupil_sampling = [n_px]);
    src.through(&mut gmt).xpupil();
    let mut pssn = CEO::<PSSN>::new().build::<ATE>(&mut src);
    let mut fwhm = ceo::Fwhm::new();
    fwhm.build(&mut src);
    fwhm.upper_bracket = 2f64 / pssn.r0() as f64;
    let mut atm = ceo::Atmosphere::new();
    let secz = 1f32 / 30f32.to_radians().cos();
    let turb_cn2_height = [275, 425, 1250, 4000, 8000, 13000]
        .iter()
        .map(|x| *x as f32 * secz)
        .collect::<Vec<f32>>();
    let turb_cn2_xi0 = [0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751];
    let wind_speed = [5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187]
        .iter()
        .map(|x| *x as f32 / secz)
        .collect::<Vec<f32>>();
    let wind_direction = [0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462];
    atm.raytrace_build(
        pssn.r0(),
        pssn.oscale,
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
    let cfd_case = "b2019_0z_0az_os_7ms";
    let duration = 200;
    let rate = 40;
    let mut ds = cfd::DomeSeeing::new(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020",
        &cfd_case,
        duration,
        Some(rate),
    );
    ds.get_keys().await?.load_opd().await?;

    //    let tau = 0.1;
    let mut xe = vec![0.; gmt2wfs.n_mode];
    let mut x = vec![0.; n_kl * 7];
    let gain = 0.5;
    gmt2wfs.qr();
    let bar = ProgressBar::new(ds.n_step as u64);
    let now = Instant::now();
    loop {
        bar.inc(1);
        atm.secs = ds.current_time;

        src.through(&mut gmt)
            .xpupil()
            .through(&mut ds)
            .through(&mut atm)
            .through(&mut pssn);
        /*
        println!(
            " {:06.3} [#{:03}]: {:6.0} [{:}]",
            ds.current_time,
            k,
            src.wfe_rms_10e(-9)[0],
            src.segment_wfe_rms_10e(-9)
                .iter()
                .map(|x| format!("{:5.0}", x))
                .collect::<Vec<String>>()
                .as_slice()
                .join(",")
        );
        */

        wfs.reset();
        gs.through(&mut gmt)
            .xpupil()
            .through(&mut ds)
            .through(&mut atm)
            .through(&mut wfs);
        wfs.process();
        xe = gmt2wfs.solve(&mut wfs.get_data());
        (0..7).for_each(|k| xe.insert(k * n_kl, 0.));
        x.iter_mut().zip(xe.iter()).for_each(|t| {
            *t.0 -= gain * *t.1 as f64;
        });
        gmt.set_m2_modes(&mut x);
        match ds.next() {
            Some(p) => {
                if p == rate {
                    //println!("Step #{}: reset PSSn!", p);
                    pssn.reset();
                }
            }
            None => break,
        };
    }
    bar.finish();
    println!(
        "Elapsed time: {:.3}s",
        now.elapsed().as_millis() as f64 * 1e-3
    );
    println!("PSSN: {:.5}", pssn.peek().estimates[0]);
    let src_fwhm = fwhm.from_complex_otf(&pssn.telescope_error_otf());
    println!("FWHM: {:.3}", src_fwhm[0].to_arcsec());
    Ok(())
}
