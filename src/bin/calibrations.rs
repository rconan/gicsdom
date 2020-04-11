use gicsdom::agws::probe::{Sensor, SH48};
use gicsdom::ceo;
use ndarray::{Array, Array2, ShapeBuilder};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::time::Instant;

fn main() {
    let mut sh48 = SH48::new();
    sh48.build(vec![0f32], vec![0f32], vec![0f32]);
    sh48.through();
    println!("WFE RMS: {:?}nm", sh48.guide_star.wfe_rms_10e(-9));

    let mut cog = ceo::Centroiding::new();
    cog.build(sh48.optics.n_side_lenslet as u32,None);
    cog.process(&sh48.sensor, None);
    println!("# valid lenslets: {}", cog.set_valid_lenslets(Some(0.9),None));

    print!("Calibrating ");
    let m1_n_mode = 17;
    let now = Instant::now();
    let mut m2_tt = ceo::Calibration::new(sh48.optics,None);
    let calibration = m2_tt.build(&cog.valid_lenslets, Some(m1_n_mode), None).calibrate(
        vec![ceo::calibrations::Mirror::M2],
        vec![
            ceo::calibrations::RigidBodyMotion::Rxyz(0, 1e-6),
            ceo::calibrations::RigidBodyMotion::Rxyz(1, 1e-6),
        ],
    );
    println!(" in {}s", now.elapsed().as_millis());
    println!(" calibration sum: {}", calibration.iter().sum::<f32>());
    let m = ceo::calibrations::pseudo_inverse(calibration, m2_tt.n_data, m2_tt.n_mode);

    /*
    let n = 2 * cog.n_valid_lenslet as usize;
    let mut d = Array2::from_shape_vec((n, 14).strides((1, n)), calibration).unwrap();
    let (u, sig, v_t) = d.svddc_inplace(UVTFlag::Some).unwrap();
    println!("eigen values: {}", sig);

    let i_sig = sig.mapv(|x| 1.0 / x);

    let l_sv = Array2::from_diag(&i_sig);
    print!("Computing the pseudo-inverse");
    let now = Instant::now();
    let m: Array2<f32> = v_t.unwrap().t().dot(&l_sv.dot(&u.unwrap().t()));
    println!(" in {}ms", now.elapsed().as_millis());
    */

    let mut m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
    m2_rbm[0][3] = 1e-6;
    m2_rbm[0][4] = 1e-6;
    m2_rbm[2][3] = 1e-6;
    m2_rbm[4][4] = 1e-6;
    sh48.gmt.update(None, Some(&m2_rbm));
    sh48.guide_star
        .through(&mut sh48.gmt)
        .xpupil()
        .lenslet_gradients(sh48.optics.n_side_lenslet, sh48.optics.lenslet_size, &mut cog);
    let centroids = cog.grab().valids(None);
    let slopes = Array::from_shape_vec((centroids.len(), 1), centroids).unwrap();
    let m2_tt = m.dot(&slopes) * 1e6;
    println!("M2 TT:\n{:+3.2}", m2_tt.into_shape((7, 2)).unwrap());
}
