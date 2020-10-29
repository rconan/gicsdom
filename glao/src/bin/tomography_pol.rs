use ceo::{
    calibrations::Mirror, calibrations::Segment, ceo, cu::Single, element::*,
    pssn::AtmosphereTelescopeError as PE, shackhartmann::Geometric, CEOInit, Calibration,
    Conversion, Cu, CEO,
};
use cfd;
use log::LevelFilter;
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use simple_logger::SimpleLogger;
use std::{env, fs::File, time::Instant};

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Results {
    time: Vec<f64>,
    wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>,
    atm_fwhm_x: f64,
    fwhm: Vec<f64>,
    pssn: Vec<f32>,
}

enum Loop {
    Closed,
    Open,
}
enum Estimator {
    LMMSE,
    LSQ,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .with_module_level("cirrus", LevelFilter::Info)
        .with_module_level("glao", LevelFilter::Info)
        .with_module_level("ceo", LevelFilter::Trace)
        .init()
        .unwrap();

    let job_idx = env::var("AWS_BATCH_JOB_ARRAY_INDEX")
        .expect("AWS_BATCH_JOB_ARRAY_INDEX env var missing")
        .parse::<usize>()
        .expect("AWS_BATCH_JOB_ARRAY_INDEX parsing failed");
    /*let index_offset = env::var("N_INDEX_OFFSET").unwrap_or("N_0".to_owned())[2..]
    .parse::<usize>()
    .expect("N_INDEX_OFFSET parsing failed");*/
    let glao_loop = match env::var("LOOP")
        .expect("LOOP:[OPEN,CLOSED] env var missing")
        .as_str()
    {
        "OPEN" => Some(Loop::Open),
        "CLOSED" => Some(Loop::Closed),
        _ => None,
    }
    .expect("LOOP is either OPEN or CLOSED");

    let duration = env::var("DURATION")
        .expect("DURATION env var missing")
        .parse::<usize>()
        .expect("DURATION parsing failed!");

    //let _science_field = env::var("SCIENCE").expect("SCIENCE:[ONAXIS,DELAUNAY21] env var mission");
    let rate = 54;
    let upload_results = true;

    let cfd_case = &cfd::get_cases()?[job_idx];
    println!("CFD CASE: {} with {} duration", cfd_case, duration);
    let mut ds = cfd::DomeSeeing::new(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020",
        &cfd_case,
        duration,
        Some(rate),
    );
    ds.get_keys().await?.load_opd().await?;

    let n_actuator = 49;
    let n_kl = 70;

    let atm_blueprint = CEO::<ATMOSPHERE>::new();
    let wfs_blueprint = CEO::<SH48>::new();
    let gs_blueprint = wfs_blueprint.guide_stars().set_on_ring(6f32.from_arcmin());
    let src_blueprint = CEO::<SOURCE>::new().set_pupil_sampling(n_actuator);

    let gmt_blueprint = CEO::<GMT>::new().set_m2_n_mode(n_kl);

    let mut lmmse = CEO::<LMMSE>::new()
        .set_atmosphere(&atm_blueprint)
        .set_guide_star(&gs_blueprint)
        .set_mmse_star(&src_blueprint)
        .set_n_side_lenslet(n_actuator - 1)
        .set_fov_diameter(10f64.from_arcmin())
        .build();
    let kln = lmmse.calibrate_karhunen_loeve(n_kl, None, None);

    let mut gmt = gmt_blueprint.build();
    let mut gs = gs_blueprint.build();
    let mut wfs = wfs_blueprint.build::<Geometric>();

    let wfs_intensity_threshold = 0.5;
    gs.through(&mut gmt).xpupil();
    wfs.calibrate(&mut gs, wfs_intensity_threshold);

    let n_px = 769;
    let mut src = ceo!(FIELDDELAUNAY21, set_pupil_sampling = [n_px]);
    src.through(&mut gmt).xpupil();

    let mut pssn = CEO::<PSSN>::new().build::<PE>(&mut src);
    let mut fwhm = ceo::Fwhm::new();
    fwhm.build(&mut src);
    fwhm.upper_bracket = 2f64 / pssn.r0() as f64;
    let atm_fwhm_x =
        ceo::Fwhm::atmosphere(500e-9, pssn.r0() as f64, pssn.oscale as f64).to_arcsec();
    let atm_fwhm_n = fwhm.from_complex_otf(&pssn.atmosphere_otf());
    println!(
        "Atm. FWHM [arcsec]: {:.3}/{:.3}",
        atm_fwhm_x,
        atm_fwhm_n[0].to_arcsec()
    );

    let mut atm = atm_blueprint
        .remove_turbulence_layer(0)
        //.set_single_turbulence_layer(0f32, Some(7f32), Some(0f32))
        .set_ray_tracing(
            25.5,
            n_px as i32,
            20f32.from_arcmin(),
            20f32,
            Some("glao_fiducial_atmosphere.bin".to_owned()),
            Some(20),
        )
        .build();

    let closed_loop_estimator = Estimator::LMMSE;
    lmmse.set_n_iteration(5);
    let wfs_rate = 18;
    let rate_ratio = rate / wfs_rate;
    let gain = 0.5;
    let mut kl_coefs = vec![0f64; n_kl * 7];
    let mut t: Vec<f64> = vec![];
    let mut wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = vec![];
    let now = Instant::now();
    match glao_loop {
        Loop::Open => {
            println!("Running open loop!");
            loop {
                atm.secs = ds.current_time;

                src.through(&mut gmt)
                    .xpupil()
                    .through(&mut ds)
                    .through(&mut atm)
                    .through(&mut pssn);
                //println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
                //println!("WFE RMS: {:?}nm", src.segment_wfe_rms_10e(-9));

                match ds.next() {
                    Some(p) => {
                        if p == rate {
                            println!("Step #{}: reset PSSn!", p);
                            //pssn.reset();
                        }
                    }
                    None => break,
                };
            }
        }
        Loop::Closed => {
            println!("Running closed loop ...");
            match closed_loop_estimator {
                Estimator::LSQ => {
                    println!("... with LSQ estimator");
                    let now = Instant::now();
                    let mut m2_2_wfs =
                        Calibration::new(&gmt_blueprint, &gs_blueprint, &wfs_blueprint);
                    m2_2_wfs.calibrate(
                        vec![Mirror::M2MODES],
                        vec![vec![Segment::Modes(1e-6, 1..n_kl)]; 7],
                        Some(wfs_intensity_threshold),
                    );
                    m2_2_wfs.qr();
                    println!("M2 calibration in {}s", now.elapsed().as_secs());
                    loop {
                        atm.secs = ds.current_time;

                        src.through(&mut gmt)
                            .xpupil()
                            .through(&mut ds)
                            .through(&mut atm)
                            .through(&mut pssn);
                        //println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
                        //println!("WFE RMS: {:?}nm", src.segment_wfe_rms_10e(-9));

                        for _ in 0..rate_ratio {
                            wfs.reset();
                            gs.through(&mut gmt)
                                .xpupil()
                                .through(&mut ds)
                                .through(&mut atm)
                                .through(&mut wfs);
                            match ds.next() {
                                Some(p) => {
                                    if p == rate {
                                        println!("Step #{}: reset PSSn!", p);
                                        //pssn.reset();
                                    }
                                }
                                None => break,
                            };
                        }
                        wfs.process();

                        let mut kl_coefs_residuals: Vec<f32> =
                            m2_2_wfs.solve(&mut wfs.get_data()).into();
                        (0..7).for_each(|k| kl_coefs_residuals.insert(k * n_kl, 0f32));
                        kl_coefs
                            .iter_mut()
                            .zip(kl_coefs_residuals.into_iter())
                            .for_each(|a| {
                                *a.0 = -gain * a.1 as f64;
                            });
                        gmt.set_m2_modes(&mut kl_coefs);
                    }
                }
                Estimator::LMMSE => {
                    println!("... with LMMSE estimator");
                    let now = Instant::now();
                    let mut m2_2_wfs =
                        Calibration::new(&gmt_blueprint, &gs_blueprint, &wfs_blueprint);
                    m2_2_wfs.calibrate(
                        vec![Mirror::M2MODES],
                        vec![vec![Segment::Modes(1e-6, 0..n_kl)]; 7],
                        None,
                    );
                    println!("M2 calibration in {}s", now.elapsed().as_secs());
                    loop {
                        t.push(ds.current_time);
                        atm.secs = ds.current_time;

                        let wfe_rms = src
                            .through(&mut gmt)
                            .xpupil()
                            .through(&mut ds)
                            .through(&mut atm)
                            .through(&mut pssn)
                            .wfe_rms_10e(-9);
                        wfe.push((
                            wfe_rms,
                            src.segment_wfe_rms_10e(-9),
                            src.segment_piston_10e(-9),
                        ));

                        let mut frame = rate_ratio;
                        let state = loop {
                            wfs.reset();
                            gs.through(&mut gmt)
                                .xpupil()
                                .through(&mut ds)
                                .through(&mut atm)
                                .through(&mut wfs);
                            match ds.next() {
                                None => break None,
                                Some(p) => {
                                    frame -= 1;
                                    if frame == 0 {
                                        break Some(p);
                                    }
                                }
                            }
                        };
                        match state {
                            Some(p) => {
                                if p == rate {
                                    println!("Step #{}: reset PSSn!", p);
                                    //pssn.reset();
                                }
                            }
                            None => break,
                        };
                        wfs.process();

                        wfs.centroids -= &m2_2_wfs.poke
                            * &Cu::<Single>::from(
                                kl_coefs.iter().map(|x| *x as f32).collect::<Vec<f32>>(),
                            );
                        lmmse.get_wavefront_estimate(&mut wfs);
                        let kl_coefs_residuals = lmmse.get_karhunen_loeve_coefficients(&kln, None);

                        kl_coefs
                            .iter_mut()
                            .zip(kl_coefs_residuals.iter())
                            .for_each(|a| {
                                *a.0 = (1f64 - gain) * *a.0 - gain * a.1;
                            });
                        gmt.set_m2_modes(&mut kl_coefs);
                    }
                }
            }
        }
    }
    println!(
        " [{}steps/{:.3}s] in {}s",
        ds.n_step,
        ds.current_time,
        now.elapsed().as_secs()
    );
    if upload_results {
        let atm_fwhm_x =
            ceo::Fwhm::atmosphere(pssn.wavelength as f64, pssn.r0() as f64, pssn.oscale as f64)
                .to_arcsec();
        let glao_fwhm = fwhm.from_complex_otf(&pssn.telescope_error_otf());
        let results = Results {
            time: t,
            wfe,
            atm_fwhm_x,
            fwhm: glao_fwhm
                .into_iter()
                .map(|x| x.to_arcsec())
                .collect::<Vec<f64>>(),
            pssn: pssn.peek().estimates.clone(),
        };
        let key = format!(
            "{}/{}/glao_{}_{}.pkl",
            "Baseline2020",
            cfd_case,
            match glao_loop {
                Loop::Open => "open",
                Loop::Closed => "closed",
            },
            match closed_loop_estimator {
                Estimator::LMMSE => "lmmse",
                Estimator::LSQ => "lsq",
            }
        );
        println!("Uploading results to: s3://gmto.modeling/{}", key);
        cirrus::dump("us-west-2", "gmto.modeling", &key, &results).await?;
    } else {
        println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
        println!("WFE RMS: {:?}nm", src.segment_wfe_rms_10e(-9));
        pssn.peek();
        println!("PSSn: {}", pssn);
        let tomo_fwhm = fwhm.from_complex_otf(&pssn.telescope_error_otf());
        println!(
            "FWHM: {:?}",
            tomo_fwhm
                .iter()
                .map(|x| x.to_arcsec())
                .collect::<Vec<f64>>()
        );
        let mut file = File::create("phase.pkl").unwrap();
        pickle::to_writer(&mut file, &src.phase(), true).unwrap();
        println!("OPEN LOOP STATE:");
        gmt.reset();
        src.through(&mut gmt)
            .xpupil()
            .through(&mut ds)
            .through(&mut atm);
        println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
        println!("WFE RMS: {:?}nm", src.segment_wfe_rms_10e(-9));
    }
    Ok(())
}
