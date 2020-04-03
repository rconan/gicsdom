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
    cog.build(sh48.n_side_lenslet as u32);
    cog.process(&sh48.sensor, None);
    println!("# valid lenslets: {}", cog.set_valid_lenslets(0.9));
    println!("c sum: {}", cog.grab().centroids.iter().sum::<f32>());

    sh48.guide_star
        .lenslet_gradients(sh48.n_side_lenslet, sh48.lenslet_size, &mut cog);
    //let c0 = cog.centroids();
    //let vc0 = cog.valid_centroids();
    println!("c sum : {}", cog.grab().centroids.iter().sum::<f32>());
    //    println!("vc sum: {}",cog.valid_centroids().iter().sum::<f32>());
    //    println!("c sum : {}",c0.iter().sum::<f32>());
    //    println!("vc sum: {}",vc0.iter().sum::<f32>());

    /*
    let mut m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];

    m2_rbm[0][3] = 1e-6;
    sh48.gmt.update(None, Some(&m2_rbm));

    sh48.guide_star.through(&mut sh48.gmt).xpupil();
    sh48.guide_star
        .lenslet_gradients(sh48.n_side_lenslet, sh48.lenslet_size, &mut cog);
    let c1 = cog.grab().valids(None);
    println!("c sum : {}", c1.iter().sum::<f32>());

    m2_rbm[0][3] = -1e-6;
    sh48.gmt.update(None, Some(&m2_rbm));
    sh48.guide_star.through(&mut sh48.gmt).xpupil();
    sh48.guide_star
        .lenslet_gradients(sh48.n_side_lenslet, sh48.lenslet_size, &mut cog);
    let c2 = cog.grab().valids(None);
    println!("c sum : {}", c2.iter().sum::<f32>());

    println!(
        "diff: {}",
        c1.iter()
            .zip(c2.iter())
            .map(|x| 0.5 * (x.0 - x.1) * 1e6)
            .sum::<f32>()
    );
     */

    let mut m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];

    let now = Instant::now();
    let mut calibration: Vec<f32> =
        Vec::with_capacity(7 * 2 * 2 * cog.n_valid_lenslet as usize * 2);
    print!("Calibrating ");
    for k in 0..7 {
        for idx in 3..5 {
            print!(".");
            // PUSH
            m2_rbm[k][idx] = 1e-6;
            sh48.gmt.update(None, Some(&m2_rbm));
            sh48.guide_star
                .through(&mut sh48.gmt)
                .xpupil()
                .lenslet_gradients(sh48.n_side_lenslet, sh48.lenslet_size, &mut cog);
            let c_push = cog.grab().valids(None);
            // PULL
            m2_rbm[k][idx] = -1e-6;
            sh48.gmt.update(None, Some(&m2_rbm));
            sh48.guide_star
                .through(&mut sh48.gmt)
                .xpupil()
                .lenslet_gradients(sh48.n_side_lenslet, sh48.lenslet_size, &mut cog);
            let c_pull = cog.grab().valids(None);
            // DIFF
            calibration.extend::<Vec<f32>>(
                c_push
                    .iter()
                    .zip(c_pull)
                    .map(|x| 0.5 * (x.0 - x.1) * 1e6)
                    .collect(),
            );
            // RESET
            m2_rbm[k][idx] = 0.0;
        }
    }
    println!(" in {}s", now.elapsed().as_millis());
    println!(" calibration sum: {}", calibration.iter().sum::<f32>());

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

    m2_rbm[0][3] = 1e-6;
    m2_rbm[0][4] = 1e-6;
    m2_rbm[2][3] = 1e-6;
    m2_rbm[4][4] = 1e-6;
    sh48.gmt.update(None, Some(&m2_rbm));
    sh48.guide_star
        .through(&mut sh48.gmt)
        .xpupil()
        .lenslet_gradients(sh48.n_side_lenslet, sh48.lenslet_size, &mut cog);
    let centroids = cog.grab().valids(None);
    let slopes = Array::from_shape_vec((centroids.len(), 1), centroids).unwrap();
    let m2_tt = m.dot(&slopes) * 1e6;
    println!("M2 TT:\n{:+3.2}", m2_tt.into_shape((7, 2)).unwrap());
}
