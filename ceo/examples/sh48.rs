use ceo::{
    calibrations, ceo, shackhartmann::Geometric as WFS_TYPE, Builder, Calibration,
    GMT, SH48,
};
use serde_pickle as pickle;
use std::fs::File;
use std::time::Instant;

fn main() {
    let mut gmt = ceo!(GMT);
    println!("M1 mode: {}",gmt.get_m1_mode_type());
    println!("M2 mode: {}",gmt.get_m2_mode_type());
    let wfs_blueprint = SH48::<WFS_TYPE>::new().set_n_sensor(1);
    let mut gs = wfs_blueprint.guide_stars().build();
    println!("GS band: {}",gs.get_photometric_band());

    let mut gmt2wfs = Calibration::new(
        &gmt,
        &gs,
        wfs_blueprint.clone(),
    );
    let mirror = vec![calibrations::Mirror::M2];
    let segments = vec![vec![calibrations::Segment::Rxyz(1e-6, Some(0..2))]; 7];
    let now = Instant::now();
    gmt2wfs.calibrate(mirror, segments, Some(0.8));
    println!(
        "GTM 2 WFS calibration [{}x{}] in {}s",
        gmt2wfs.n_data,
        gmt2wfs.n_mode,
        now.elapsed().as_secs()
    );
    let poke_sum = gmt2wfs.poke.from_dev().iter().sum::<f32>();
    println!("Poke sum: {}", poke_sum);

    let mut file = File::create("poke.pkl").unwrap();
    pickle::to_writer(&mut file, &gmt2wfs.poke.from_dev(), true).unwrap();

    let mut wfs = wfs_blueprint.build();
    let mut src = ceo!(SOURCE);
    let mut atm = ceo!(ATMOSPHERE);


    gs.through(&mut gmt).xpupil();
    println!("GS WFE RMS: {}nm", gs.wfe_rms_10e(-9)[0]);
    wfs.calibrate(&mut gs, 0.8);
    println!("# valid lenslet: {}", wfs.n_valid_lenslet());

    gmt.set_m2_segment_state(2, &[0., 0.0, 0.], &[1e-6, 0.0, 0.]);
    gmt.set_m2_segment_state(5, &[0., 0.0, 0.], &[0., 1e-6, 0.]);
    gmt.set_m2_segment_state(7, &[0., 0.0, 0.], &[1e-6, 1e-6, 0.]);
    wfs.reset();
    gs.through(&mut gmt).xpupil().through(&mut wfs);
    wfs.process();

    let a = gmt2wfs.qr().solve(&mut wfs.get_data());
    let s: Vec<f32> = (&gmt2wfs.poke * &a).into();

    let mut file = File::create("slopes.pkl").unwrap();
    pickle::to_writer(&mut file, &(wfs.get_data().from_dev(), s), true).unwrap();

    Vec::<f32>::from(a)
        .into_iter()
        .map(|x| x * 1e6)
        .collect::<Vec<f32>>()
        .chunks(2)
        .enumerate()
        .for_each(|x| println!("#{}: [{:+0.1},{:+0.1}]", 1 + x.0, x.1[0], x.1[1]));
    //    println!("M2 TT: {:#?}", a);

    println!(
        "WFE RMS [nm] without and with atmosphere: {:.0}/{:.0}",
        src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0],
        src.through(&mut gmt)
            .xpupil()
            .through(&mut atm)
            .wfe_rms_10e(-9)[0]
    );
}
