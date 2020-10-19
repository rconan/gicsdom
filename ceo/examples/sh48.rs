use ceo::{
    calibrations, ceo, element::*, shackhartmann::Geometric as WFS_TYPE, CEOBluePrint, Calibration,
    Conversion, CEO,
};
use std::time::Instant;

fn main() {
    let gmt_blueprint = CEO::<GMT>::new();
    let wfs_blueprint = CEO::<SH48>::new().set_n_sensor(1);
    let gs_blueprint = wfs_blueprint.guide_stars();

    let mut gmt2wfs = Calibration::new(gmt_blueprint.clone(), gs_blueprint.clone());
    let mirror = vec![calibrations::Mirror::M2];
    let segments = vec![vec![calibrations::Segment::Rxyz(1e-6, Some(0..2))]; 7];
    let now = Instant::now();
    gmt2wfs.calibrate(mirror, segments, wfs_blueprint.clone());
    println!(
        "GTM 2 WFS calibration [{}x{}] in {}s",
        gmt2wfs.n_data,
        gmt2wfs.n_mode,
        now.elapsed().as_secs()
    );
    let poke_sum = gmt2wfs.poke.from_dev().iter().sum::<f32>();
    println!("Poke sum: {}",poke_sum);

    let mut gmt = gmt_blueprint.build();
    let mut wfs = wfs_blueprint.build::<WFS_TYPE>();
    let mut gs = gs_blueprint.build();
    let mut src = ceo!(SOURCE);
    let mut atm = ceo!(ATMOSPHERE);

    gs.through(&mut gmt).xpupil();
    println!("GS WFE RMS: {}nm", gs.wfe_rms_10e(-9)[0]);
    wfs.calibrate(&mut gs, 0.);
    println!("# valid lenslet: {}", wfs.n_valid_lenslet());

    gmt.set_m1_segment_state(1, &[0., 0.0, 0.], &[1e-6, 0.0, 0.]);
    wfs.reset();
    gs.through(&mut gmt).xpupil().through(&mut wfs);
    wfs.process();

    let a = gmt2wfs
        .solve(&mut wfs.get_data())
        .iter()
        .map(|x| x * 1e6)
        .collect::<Vec<f32>>();
    println!("M2 TT: {:#?}", a);

    println!(
        "WFE RSM [nm] without and with atmosphere: {:.0}/{:.0}",
        src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0],
        src.through(&mut gmt)
            .xpupil()
            .through(&mut atm)
            .wfe_rms_10e(-9)[0]
    );
}
