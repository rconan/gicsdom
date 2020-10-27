use ceo::{
    calibrations::Mirror, calibrations::Segment, ceo, cu::Single, element::*,
    shackhartmann::Geometric, CEOInit, Calibration, Conversion, Cu, CEO,
};
use serde_pickle as pickle;
use simple_logger;
use std::fs::File;

fn main() {
    simple_logger::init().unwrap();

    let n_actuator = 49;
    let n_kl = 70;

    let atm_blueprint = CEO::<ATMOSPHERE>::new()
        .set_r0_at_zenith(0.15)
        .set_oscale(60.)
        .set_zenith_angle(0.)
        .set_turbulence_profile(TurbulenceProfile {
            n_layer: 1,
            altitude: vec![0f32],
            xi0: vec![1f32],
            wind_speed: vec![0f32],
            wind_direction: vec![0f32],
        });
    let wfs_blueprint = CEO::<SH48>::new();
    let gs_blueprint = wfs_blueprint.guide_stars().set_on_ring(6f32.from_arcmin());
    let src_blueprint = CEO::<SOURCE>::new().set_pupil_sampling(n_actuator);

    let gmt_blueprint = CEO::<GMT>::new().set_m2_n_mode(n_kl);
    let mut mmse_src = src_blueprint.build();

    let mut lmmse = CEO::<LMMSE>::new()
        .set_atmosphere(&atm_blueprint)
        .set_guide_star(&gs_blueprint)
        .set_mmse_star(&src_blueprint)
        .set_n_side_lenslet(n_actuator - 1)
        .build();

    let mut gmt = gmt_blueprint.build();
    let mut atm = atm_blueprint.build();
    let mut gs = gs_blueprint.build();
    let mut wfs = wfs_blueprint.build::<Geometric>();
    let kln = lmmse.calibrate_karhunen_loeve(n_kl, None, None);

    let wfs_intensity_threshold = 0.5;
    gs.through(&mut gmt).xpupil();
    wfs.calibrate(&mut gs, wfs_intensity_threshold);

    let mut m2_2_wfs = Calibration::new(&gmt_blueprint, &gs_blueprint, &wfs_blueprint);
    m2_2_wfs.calibrate(
        vec![Mirror::M2MODES],
        vec![vec![Segment::Modes(1e-6, 0..n_kl)]; 7],
        None,
    );

    let mut src = ceo!(SOURCE);
    let mut kl_coefs = vec![0f64; n_kl * 7];

    src.through(&mut gmt).xpupil().through(&mut atm);
    println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
    println!("WFE RMS: {:?}nm", src.segment_wfe_rms_10e(-9));

    let mut data = vec![];

    for _ in 0..20 {
        wfs.reset();
        gs.through(&mut gmt)
            .xpupil()
            .through(&mut atm)
            .through(&mut wfs);
        wfs.process();
        let s_cl: Vec<f32> = wfs.centroids.from_dev();

        let s_kl = m2_2_wfs.poke.clone()
            * Cu::<Single>::from(kl_coefs.iter().map(|x| *x as f32).collect::<Vec<f32>>());
        wfs.centroids -= s_kl.clone();
        let s_ol: Vec<f32> = wfs.centroids.from_dev();
        data.push((s_cl, Vec::<f32>::from(s_kl), s_ol));
        //let mut lmmse_phase = lmmse.get_wavefront_estimate(&mut wfs).phase_as_ptr();
        lmmse.get_wavefront_estimate(&mut wfs);
        let kl_coefs_residuals = lmmse.get_karhunen_loeve_coefficients(&kln, None);

        let gain = 0.5;
        kl_coefs
            .iter_mut()
            .zip(kl_coefs_residuals.iter())
            .for_each(|a| {
                *a.0 = (1f64 - gain) * *a.0 - gain * a.1;
            });
        gmt.set_m2_modes(&mut kl_coefs);

        println!("# of iteration: {}", lmmse.get_n_iteration());
        src.through(&mut gmt).xpupil().through(&mut atm);
        println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
        println!("WFE RMS: {:?}nm", src.segment_wfe_rms_10e(-9));
        //let src_phase = mmse_src.phase().clone();
        //mmse_src.sub(&mut lmmse_phase);
        //println!("Residual WFE RMS: {}nm", mmse_src.wfe_rms_10e(-9)[0]);
    }

    let mut file = File::create("slopes.pkl").unwrap();
    pickle::to_writer(&mut file, &data, true).unwrap();

    /*
    let phase = Vec::<f32>::from(lmmse_phase);
    let mut file = File::create("tomography.pkl").unwrap();
    pickle::to_writer(&mut file, &(src_phase, phase), true).unwrap();

    let src_phase: Vec<f32> = src.phase_as_ptr().into();
    let mut file = File::create("SRC_phase.pkl").unwrap();
    pickle::to_writer(&mut file, &src_phase, true).unwrap();

    src.through(&mut gmt).xpupil().through(&mut atm);
    println!("KL residual WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);

    let kl_phase: Vec<f32> = src.phase_as_ptr().into();
    let mut file = File::create("KL_phase.pkl").unwrap();
    pickle::to_writer(&mut file, &kl_phase, true).unwrap();
    */
}
