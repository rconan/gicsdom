use gicsdom::agws::probe::{GmtState, Message};
use gicsdom::agws::probe::{Sensor, SH48};
use gicsdom::ceo;
use gicsdom::ceo::Conversion;
use ndarray;
use ndarray::{Array, Array2, ShapeBuilder};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::ops::Range;
use std::time::Instant;
use gicsdom::astrotools;
use gicsdom::agws;

fn main() {
    let mut sh48 = SH48::new();
    sh48.build(vec![0f32], vec![0f32], vec![0f32]);
    sh48.through(None);
    println!("WFE RMS: {:?}nm", sh48.guide_star.wfe_rms_10e(-9));

    let mut cog = ceo::Centroiding::new();
    cog.build(sh48.optics.n_side_lenslet as u32, None);
    cog.process(&sh48.sensor, None);
    println!(
        "# valid lenslets: {}",
        cog.set_valid_lenslets(Some(0.9), None)
    );

    let mut cog1 = ceo::Centroiding::new();
    cog1.build(sh48.optics.n_side_lenslet as u32, None);
    sh48.through(None);
    cog1.process(&sh48.sensor, Some(&cog));
    let centroids = cog1.grab().valids(Some(&cog.valid_lenslets));
    println!("centroids sum: {}", centroids.iter().sum::<f32>());

    print!("Calibrating ");
    let m1_n_mode = 17;
    let now = Instant::now();
    let mut rbm_calib = ceo::Calibration::new(sh48.optics, None);
    let mut bm_calib = ceo::Calibration::new(sh48.optics, None);

    let mirrors = vec![ceo::calibrations::Mirror::M2];
    let rbms = (0..7)
        .into_iter()
        .map(|k| {
            if k == 6 {
                vec![
                    //  ceo::calibrations::Segment::Txyz(1e-6, None),
                    ceo::calibrations::Segment::Rxyz(1e-6, Some(0..2)),
                ]
            } else {
                vec![
                    // ceo::calibrations::Segment::Txyz(1e-6, None),
                    ceo::calibrations::Segment::Rxyz(1e-6, Some(0..2)),
                ]
            }
        })
        .collect::<Vec<Vec<ceo::calibrations::Segment>>>();

    //let mirrors = vec![ceo::calibrations::Mirror::M1MODES];
    let modes = vec![vec![ceo::calibrations::Segment::Modes(1e-6, 0..17)]; 7];

    let rbm_calibration = rbm_calib
        .build(0., 0., &cog.valid_lenslets, Some(m1_n_mode), None)
        .calibrate(mirrors.clone(), rbms.clone());
    let bm_calibration = bm_calib
        .build(0., 0., &cog.valid_lenslets, Some(m1_n_mode), None)
        .calibrate(vec![ceo::calibrations::Mirror::M1MODES], modes.clone());
    println!(" in {}s", now.elapsed().as_millis());
    println!(
        "Calibration RBM matrix: [{};{}]",
        rbm_calib.n_data, rbm_calib.n_mode
    );
    println!(
        "Calibration BM matrix: [{};{}]",
        bm_calib.n_data, bm_calib.n_mode
    );
    //let calibration = vec![rbm_calibration, bm_calibration].into_iter().flatten().collect::<Vec<f32>>();
    //let n_mode = vec![rbm_calib.n_mode,bm_calib.n_mode].iter().sum::<u32>();
    //println!(" calibration sum: {}", calibration.iter().sum::<f32>());
    let mut d = ceo::calibrations::composite(
        vec![rbm_calibration, bm_calibration],
        vec![rbm_calib.n_mode, bm_calib.n_mode],
        Some(1)
    );
    let m = ceo::calibrations::pseudo_inverse(
        &mut d,
        None
    );

    let mut m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
    let a = 1e-6;
    println!("a:{:.0}mas", a.to_mas());


    m2_rbm[0][3] = a;
    m2_rbm[0][4] = a;
    m2_rbm[2][3] = a;
    m2_rbm[4][4] = a;
    m2_rbm[6][3] = a;
    m2_rbm[6][4] = a;

    //    sh48.gmt.update(None, Some(&m2_rbm), None);

    let mut m1_mode = vec![vec![0f64; 17]; 7];
    m1_mode[0][0] = 100e-9;
    m1_mode[2][8] = -30e-9;
    m1_mode[4][9] = 60e-9;
    sh48.gmt.update(None, Some(&m2_rbm), Some(&m1_mode));

    sh48.guide_star
        .through(&mut sh48.gmt)
        .xpupil()
        .lenslet_gradients(
            sh48.optics.n_side_lenslet,
            sh48.optics.lenslet_size,
            &mut cog,
        );
    let centroids = cog.grab().valids(None);
    println!("centroids sum: {}", centroids.iter().sum::<f32>());
    let slopes = Array::from_shape_vec((centroids.len(), 1), centroids).unwrap();
    let coms = m.dot(&slopes);
    println! {"{}",coms.sum()*1e9};

    let mut state = GmtState::new(Some(17));
    state.fill_from_array(mirrors.clone(), rbms.clone(), &coms);
    let q = coms.slice(ndarray::s![14..,..]).to_owned();
    state.fill_from_array(vec![ceo::calibrations::Mirror::M1MODES], modes.clone(),
                          &q);
    println!(
        "{}",state
    );
    /*
       println!(
       "M2 RBM:\n{:.0}",
       1e6*m2_tt.into_shape((m2_calib.n_mode as usize / 2, 2)).unwrap()
       );
     */

}
