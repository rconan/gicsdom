use ceo::{ceo, element::*, shackhartmann::Geometric, CEOInit, Conversion, Mask, CEO};
use serde_pickle as pickle;
use std::fs::File;

fn main() {
    let n_actuator = 49;
    let n_kl = 70;

    let atm_blueprint = CEO::<ATMOSPHERE>::new()
        .set_r0_at_zenith(0.15)
        .set_oscale(30.)
        .set_zenith_angle(0.)
        .set_turbulence_profile(TurbulenceProfile {
            n_layer: 1,
            altitude: vec![0f32],
            xi0: vec![1f32],
            wind_speed: vec![0f32],
            wind_direction: vec![0f32],
        });
    let wfs_blueprint = CEO::<SH48>::new().set_n_sensor(3);
    let gs_blueprint = wfs_blueprint.guide_stars().set_on_ring(6f32.from_arcmin());
    let src_blueprint = CEO::<SOURCE>::new().set_pupil_sampling(n_actuator);

    let mut gmt = ceo!(GMT, set_m2_n_mode = [n_kl]);
    let mut mmse_src = src_blueprint.build();

    mmse_src.through(&mut gmt).xpupil();
    let mut dm = Mask::new();
    dm.build(n_actuator * n_actuator)
        .filter(&mut mmse_src.amplitude().into());

    let mut kl: Vec<Vec<f32>> = vec![];
    for s in 0..7 {
        for k in 1..n_kl {
            gmt.set_m2_modes_ij(s, k, 1e-6);
            mmse_src.through(&mut gmt).xpupil();
            kl.push(mmse_src.phase_as_ptr().into());
        }
    }
    gmt.reset();
    let kl_norm = kl
        .iter()
        .map(|x| x.iter().map(|y| (y * y) as f64).sum::<f64>())
        .collect::<Vec<f64>>();
    let mut file = File::create("KL.pkl").unwrap();
    pickle::to_writer(&mut file, &kl, true).unwrap();

    let mut lmmse = CEO::<LMMSE>::new()
        .set_atmosphere(&atm_blueprint)
        .set_guide_star(&gs_blueprint)
        .set_mmse_star(&src_blueprint)
        .set_n_side_lenslet(n_actuator - 1)
        .set_pupil_mask(dm)
        .build();

    let mut atm = atm_blueprint.build();
    let mut gs = gs_blueprint.build();
    let mut wfs = wfs_blueprint.build::<Geometric>();

    gs.through(&mut gmt).xpupil();
    wfs.calibrate(&mut gs, 0.5);

    wfs.reset();
    gs.through(&mut gmt)
        .xpupil()
        .through(&mut atm)
        .through(&mut wfs);
    wfs.process();
    let mut lmmse_phase = lmmse.get_wavefront_estimate(&mut wfs).phase_as_ptr();
    println!("# of iteration: {}", lmmse.get_n_iteration());
    mmse_src.through(&mut gmt).xpupil().through(&mut atm);
    println!("WFE RMS: {}nm", mmse_src.wfe_rms_10e(-9)[0]);
    let src_phase = mmse_src.phase().clone();
    mmse_src.sub(&mut lmmse_phase);
    println!("Residual WFE RMS: {}nm", mmse_src.wfe_rms_10e(-9)[0]);

    let phase = Vec::<f32>::from(lmmse_phase);
    let mut kl_coefs = kl
        .iter()
        .zip(kl_norm.iter())
        .map(|x| {
            x.0.iter()
                .zip(phase.iter())
                .map(|y| -1e-6 * (y.0 * y.1) as f64)
                .sum::<f64>()
                / x.1
        })
        .collect::<Vec<f64>>();
    (0..7).for_each(|k| kl_coefs.insert(k * n_kl, 0.));

    let mut file = File::create("KL_coefs.pkl").unwrap();
    pickle::to_writer(&mut file, &kl_coefs, true).unwrap();

    let mut file = File::create("tomography.pkl").unwrap();
    pickle::to_writer(&mut file, &(src_phase, phase), true).unwrap();

    let mut src = ceo!(SOURCE);

    src.through(&mut gmt).xpupil().through(&mut atm);
    println!("WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);
    let src_phase: Vec<f32> = src.phase_as_ptr().into();
    let mut file = File::create("SRC_phase.pkl").unwrap();
    pickle::to_writer(&mut file, &src_phase, true).unwrap();

    gmt.set_m2_modes(&mut kl_coefs);
    src.through(&mut gmt).xpupil().through(&mut atm);
    println!("KL residual WFE RMS: {}nm", src.wfe_rms_10e(-9)[0]);

    let kl_phase: Vec<f32> = src.phase_as_ptr().into();
    let mut file = File::create("KL_phase.pkl").unwrap();
    pickle::to_writer(&mut file, &kl_phase, true).unwrap();

}
