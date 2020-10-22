use ceo::{
    calibrations, ceo, element::*, shackhartmann::Geometric as WFS_TYPE, CEOInit, Calibration,
    Conversion, CEO,
};
//use serde_pickle as pickle;
//use std::fs::File;
use std::time::Instant;

fn main() {
    let gmt_blueprint = CEO::<GMT>::new().set_m2_n_mode(70);
    let wfs_blueprint = CEO::<SH48>::new();
    let gs_blueprint = wfs_blueprint.guide_stars().set_on_ring(6f32.from_arcmin());

    let mut gmt2wfs = Calibration::new(&gmt_blueprint, &gs_blueprint, &wfs_blueprint);
    let mirror = vec![calibrations::Mirror::M2MODES];
    let segments = vec![vec![calibrations::Segment::Modes(1e-6, 1..70)]; 7];
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

    //    let mut file = File::create("poke.pkl").unwrap();
    //    pickle::to_writer(&mut file, &gmt2wfs.poke.from_dev(), true).unwrap();

    let mut gmt = gmt_blueprint.build();
    let mut wfs = wfs_blueprint.build::<WFS_TYPE>();
    let mut gs = gs_blueprint.build();
    let mut src = ceo!(SOURCE);
    let mut atm = ceo!(ATMOSPHERE);

    gs.through(&mut gmt).xpupil();
    println!("GS WFE RMS: {}nm", gs.wfe_rms_10e(-9)[0]);
    wfs.calibrate(&mut gs, 0.8);
    println!("# valid lenslet: {}", wfs.n_valid_lenslet());

    let mut m = vec![vec![0f64; 70]; 7];
    m[0][1] = 1e-6;
    m[1][3] = -1e-6;
    m[3][4] = 1e-6;
    m[6][2] = 1e-6;
    m[6][4] = -1e-6;
    gmt.update(None, None, None, Some(&m));
    wfs.reset();
    gs.through(&mut gmt).xpupil().through(&mut wfs);
    wfs.process();

    //    let mut file = File::create("slopes.pkl").unwrap();
    //    pickle::to_writer(&mut file, &wfs.get_data().from_dev(), true).unwrap();

    let a = gmt2wfs
        .qr()
        .solve(&mut wfs.get_data())
        .iter()
        .map(|x| x * 1e6)
        .collect::<Vec<f32>>()
        .chunks(69)
        .enumerate()
        .for_each(|x| {
            println!(
                "#{}: [{}]",
                1 + x.0,
                x.1.iter()
                    .take(5)
                    .map(|x| format!("{:+0.1}", x))
                    .collect::<Vec<String>>()
                    .as_slice()
                    .join(",")
            )
        });
    //    println!("M2 TT: {:#?}", a);

    /*
    println!(
        "WFE RSM [nm] without and with atmosphere: {:.0}/{:.0}",
        src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0],
        src.through(&mut gmt)
            .xpupil()
            .through(&mut atm)
            .wfe_rms_10e(-9)[0]
    );
     */
}
