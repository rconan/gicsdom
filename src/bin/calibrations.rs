use gicsdom::agws::probe::{GmtState, Message};
use gicsdom::agws::probe::{Sensor, SH48};
use gicsdom::ceo;
use gicsdom::ceo::Conversion;
use ndarray::{Array, Array2, ShapeBuilder};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::ops::Range;
use std::time::Instant;

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

    print!("Calibrating ");
    let m1_n_mode = 17;
    let now = Instant::now();
    let mut m2_calib = ceo::Calibration::new(sh48.optics, None);
    let mirrors = vec![ceo::calibrations::Mirror::M1,ceo::calibrations::Mirror::M2];
    let rbms = (0..7)
        .into_iter()
        .map(|k| {
            if k == 6 {
                vec![
                  //  ceo::calibrations::RigidBodyMotion::Txyz(1e-6, None),
                    ceo::calibrations::RigidBodyMotion::Rxyz(1e-6, Some(0..2)),
                ]
            } else {
                vec![
                   // ceo::calibrations::RigidBodyMotion::Txyz(1e-6, None),
                    ceo::calibrations::RigidBodyMotion::Rxyz(1e-6, Some(0..2)),
                ]
            }
        })
        .collect::<Vec<Vec<ceo::calibrations::RigidBodyMotion>>>();
    let calibration = m2_calib
        .build(&cog.valid_lenslets, Some(m1_n_mode), None)
        .calibrate(mirrors.clone(), rbms.clone());
    println!(" in {}s", now.elapsed().as_millis());
    println!(
        "Calibration matrix: [{};{}]",
        m2_calib.n_data, m2_calib.n_mode
    );
    println!(" calibration sum: {}", calibration.iter().sum::<f32>());
    let m = ceo::calibrations::pseudo_inverse(vec![calibration], m2_calib.n_mode as usize, None);

    let mut m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
    let a = 1e-6;
    println!("a:{:.0}mas", a.to_mas());
    /*
    m2_rbm[0][3] = a;
    m2_rbm[0][4] = a;
    m2_rbm[2][3] = a;
    m2_rbm[4][4] = a;
    */
    m2_rbm[0][3] = a;
    m2_rbm[0][4] = a;
    sh48.gmt.update(None, Some(&m2_rbm));
    sh48.guide_star
        .through(&mut sh48.gmt)
        .xpupil()
        .lenslet_gradients(
            sh48.optics.n_side_lenslet,
            sh48.optics.lenslet_size,
            &mut cog,
        );
    let centroids = cog.grab().valids(None);
    let slopes = Array::from_shape_vec((centroids.len(), 1), centroids).unwrap();
    let m2_tt = m.dot(&slopes);

    let mut state = GmtState::new();
    println!("{}", state.fill_from_array(mirrors, rbms, &m2_tt));
    /*
        println!(
        "M2 RBM:\n{:.0}",
        1e6*m2_tt.into_shape((m2_calib.n_mode as usize / 2, 2)).unwrap()
        );
    */
}
