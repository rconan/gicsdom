//extern crate ndarray;
//extern crate ndarray_linalg;

use std::time::Instant;
use std::f32;
use ndarray::{s, stack, Array, Array2, Axis};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
//use rand::distributions::{Normal, Uniform};
//use rand::thread_rng;
use gicsdom::ceo; //::{Gmt, Source,Geometric_shack_hartmann,Diffractive_shack_hartmann};
use indicatif::{ProgressBar, ProgressStyle};

fn main() {
    let n_side_lenslet = 48;
    let n_px_lenslet = 16;
    //let n_px = n_side_lenslet*16 + 1;
    let pupil_size = 25.5;
    let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;
    let n_gs = 3;
    let m1_n_mode = 27;

    // GMT initialization
    let mut gmt = ceo::Gmt::new(m1_n_mode, None);
    gmt.build();

    // WFS initialization
    let d = pupil_size / n_side_lenslet as f64;
    let mut gwfs = ceo::Geometric_ShackHartmann::new(n_gs, n_side_lenslet, n_px_lenslet, d);
    gwfs.build();
    let mut wfs = ceo::Diffractive_ShackHartmann::new(n_gs, n_side_lenslet, n_px_lenslet, d);
    wfs.build(8, Some(24), None);

    // Guide stars initialization
    let z = 6. * f32::consts::PI / 180. / 60.;
    let a = 2. * f32::consts::PI / 3.;
    let mut gs = gwfs.new_guide_stars();
    gs.build(
        "V",
        vec![z, z, z],
        vec![0.0 * a, a, 2.0 * a],
        vec![0.0, 0.0, 0.0],
    );

    // Sensor calibration
    gs.through(&mut gmt);
    gwfs.calibrate(&mut gs, 0.0).unwrap();
    wfs.calibrate(&mut gs, 0.0);

    println!("WFE RMS: {:?}", gs.wfe_rms_10e(-6));

    gs.through(&mut gmt).through(&mut gwfs);
    gwfs.process();

    // Science initialization
    let mut src = ceo::Source::new(1, pupil_size, pupil_sampling);
    src.build("V", vec![0.0], vec![0.0], vec![0.0]);
    println!("WFE RMS: {:?}", src.through(&mut gmt).wfe_rms_10e(-9));

    //        println!("WFE RMS: {:?}",gs.wfe_rms_10e(-6).iter()<f32>.sum());

    //println!("{:?}",gwfs.centroids);

    // ----------------------------------------------------------------------------
    // GMT CALIBRATION
    println!("Calibrating the GMT ...");
    let now = Instant::now();
    let n_rbm: usize = 84;
    let mut c_c_p: Vec<f32> = Vec::new();
    let mut c_c_m: Vec<f32> = Vec::new();
    let n_c: usize = (n_side_lenslet.pow(2) as usize) * 2 * n_gs as usize;
    gmt.reset();
    let pb = ProgressBar::new(7);
    pb.set_style(ProgressStyle::default_bar().template("{bar:40.cyan/blue} {msg}"));
    if n_rbm > 0 {
        for mid in 1..3 {
            let m_ = String::from("M") + &mid.to_string() + " RBM";
            pb.set_message(&m_);
            pb.reset();
            for sid in 1..8 {
                pb.inc(1);
                for tr in 1..3 {
                    for a in 0..3 {
                        gmt.reset();
                        let mut t_xyz = vec![0.0; 3];
                        let mut r_xyz = vec![0.0; 3];
                        if tr == 1 {
                            t_xyz[a] = 1e-6;
                        }
                        if tr == 2 {
                            r_xyz[a] = 1e-6;
                        }
                        if mid == 1 {
                            gmt.set_m1_segment_state(sid, &t_xyz, &r_xyz);
                        }
                        if mid == 2 {
                            gmt.set_m2_segment_state(sid, &t_xyz, &r_xyz);
                        }
                        gs.through(&mut gmt).through(&mut gwfs);
                        gwfs.process();
                        c_c_p.append(&mut gwfs.centroids);
                    }
                }
            }
        }
    }

    let pb_bm = ProgressBar::new(m1_n_mode as u64);
    pb_bm.set_style(ProgressStyle::default_bar().template("{bar:40.cyan/blue} {msg}"));
    if m1_n_mode > 0 {
        //println!("Bending modes ...");
        let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
        gmt.reset();
        for sid in 0..7 {
            let s_ = String::from("S") + &(sid + 1).to_string() + "BM";
            pb_bm.set_message(&s_);
            pb_bm.reset();
            for k_a in 0..m1_n_mode {
                pb_bm.inc(1);
                let k = k_a + sid * m1_n_mode;
                a[k as usize] = 1e-6;
                gmt.set_m1_modes(&mut a);
                gs.through(&mut gmt).through(&mut gwfs);
                gwfs.process();
                c_c_p.append(&mut gwfs.centroids);
                a[k as usize] = 0.0;
                gmt.set_m1_modes(&mut a);
            }
        }
    }

    let _d_p = Array2::from_shape_vec((c_c_p.len() / n_c, n_c), c_c_p)
        .unwrap()
        .t()
        .to_owned()
        * 1e6;

    //        println!("Rigid body motion ...");
    if n_rbm > 0 {
        for mid in 1..3 {
            let m_ = String::from("M") + &mid.to_string() + " RBM";
            pb.set_message(&m_);
            pb.reset();
            for sid in 1..8 {
                pb.inc(1);
                for tr in 1..3 {
                    for a in 0..3 {
                        gmt.reset();
                        let mut t_xyz = vec![0.0; 3];
                        let mut r_xyz = vec![0.0; 3];
                        if tr == 1 {
                            t_xyz[a] = -1e-6;
                        }
                        if tr == 2 {
                            r_xyz[a] = -1e-6;
                        }
                        if mid == 1 {
                            gmt.set_m1_segment_state(sid, &t_xyz, &r_xyz);
                        }
                        if mid == 2 {
                            gmt.set_m2_segment_state(sid, &t_xyz, &r_xyz);
                        }
                        gs.through(&mut gmt).through(&mut gwfs);
                        gwfs.process();
                        //let sum2: f32 = c.iter().sum();
                        //println!("Centroids sum: {:.16}", sum2);
                        //println!("Diff. centroids sum: {:.14}", sum2 - sum1);
                        c_c_m.append(&mut gwfs.centroids);
                    }
                }
            }
        }
    }

    if m1_n_mode > 0 {
        //println!("Bending modes ...");
        let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
        gmt.reset();
        for sid in 0..7 {
            let s_ = String::from("S") + &(sid + 1).to_string() + "BM";
            pb_bm.set_message(&s_);
            pb_bm.reset();
            for k_a in 0..m1_n_mode {
                pb_bm.inc(1);
                let k = k_a + sid * m1_n_mode;
                a[k as usize] = -1e-6;
                gmt.set_m1_modes(&mut a);
                gs.through(&mut gmt).through(&mut gwfs);
                gwfs.process();
                c_c_m.append(&mut gwfs.centroids);
                a[k as usize] = 0.0;
                gmt.set_m1_modes(&mut a);
            }
        }
    }
    pb.finish();
    pb_bm.finish();

    println!("WFS centroids #     : {}", gwfs.n_centroids);
    println!("WFS centroids length: {}", gwfs.centroids.len());
    println!("CCM length: {}", c_c_m.len());

    let _d_m = Array2::from_shape_vec((c_c_m.len() / n_c, n_c), c_c_m)
        .unwrap()
        .t()
        .to_owned()
        * 1e6;
    let mut _d = (_d_p - _d_m) * 0.5;
    println!(" in {}s", now.elapsed().as_secs());

    println!("{:?}", _d);
    println!("shape={:?}, strides={:?}", _d.shape(), _d.strides());
    //        println!("d sum: {}",_d.into.sum());

    print!("SVD decomposition");
    let now = Instant::now();
    let n_sv = n_rbm + m1_n_mode as usize * 7;
    let (u, sig, v_t) = _d.svddc_inplace(UVTFlag::Some).unwrap();
    println!(" in {}ms", now.elapsed().as_millis());
    println!("Singular values:\n{}", sig);
    let mut i_sig = sig.mapv(|x| 1.0 / x);
    for k in 0..14 {
        i_sig[n_sv - k - 1] = 0.0;
    }

    let _u = u.unwrap();
    let _vt = v_t.unwrap();

    let l_sv = Array2::from_diag(&i_sig);
    let z: Array2<f32> = Array2::zeros((n_sv, n_c - n_sv));
    let ss = stack(Axis(1), &[l_sv.view(), z.view()]).unwrap();
    println!("SS shape: {:?}", ss.shape());
    print!("Computing the pseudo-inverse");
    let now = Instant::now();
    let __m = _vt.t().dot(&l_sv.dot(&_u.t()));
    println!(" in {}ms", now.elapsed().as_millis());

    // Initial conditions
    let mut t_xyz = vec![0.0; 3];
    let mut r_xyz = vec![0.0; 3];
    let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
    gmt.reset();
    t_xyz[1] = 1e-6;
    gmt.set_m1_segment_state(1, &t_xyz, &r_xyz);
    r_xyz[0] = 1e-6;
    gmt.set_m2_segment_state(2, &t_xyz, &r_xyz);
    if m1_n_mode > 0 {
        a[0] = 1e-6;
        a[5] = 1e-6;
        a[11] = 1e-6;
        gmt.set_m1_modes(&mut a);
    }
    gs.through(&mut gmt).through(&mut gwfs);
    gwfs.process();

    let slopes = Array::from_shape_vec((n_c, 1), gwfs.centroids.clone()).unwrap();
    let _c_e = __m.dot(&slopes);
    //let _c_e2 = _c_e.into_shape((14, 6)).unwrap();
    let rbm = _c_e.slice(s![..84, ..]);
    println!(
        "RBM:\n{:.2}",
        rbm.to_owned().into_shape((14, 6)).unwrap() * 1e6
    );
    if m1_n_mode > 0 {
        let bm = _c_e.slice(s![84.., ..]);
        println!(
            "BM:\n{:.2}",
            bm.to_owned().into_shape((7, m1_n_mode as usize)).unwrap() * 1e6
        );
    }

    // ------------------------------------------------------------------------------------------------------
    let closed_loop = true;
    if closed_loop {
        //        let normal: Normal = Normal::new(0,1e-7);
        //        let mut com = MatrixMN::<f64, U6, U14>::from_distribution(&normal,&mut rng);
        //let mut com: Array2<f32> = Array2::zeros((6,14));
        let gain = 0.5f32;

        let mut src_wfe_rms = 0f32;
        let mut onaxis_pssn = 0f32;
        let mut t_xyz = vec![0.0; 3];
        let mut r_xyz = vec![0.0; 3];
        let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];

        src_wfe_rms = src.through(gmt.reset()).wfe_rms_10e(-9)[0];

        let mut src_pssn = ceo::PSSn::new(15e-2, 25.0, 0.0);
        src_pssn.build(&mut src);
        println!(
            "WFE RMS: {:.3} nm; PSSn: {}",
            src_wfe_rms,
            src_pssn.reset(&mut src)
        );

        //println!("{:.2}", com[(0,0)] * 1e6);
        let mut id = 0;
        let px = (pupil_size / (pupil_sampling as f64 - 1.0)) as f32;

        // GMT initial conditions
        let mut rbm = Array2::<f32>::zeros((14, 6));
        rbm[(0, 0)] = 1e-5;
        rbm[(8, 4)] = 1e-6;
        rbm[(2, 3)] = -1e-6;
        rbm[(11, 1)] = -1e-5;
        let mut bm = Array2::<f32>::zeros((7, m1_n_mode as usize));
        bm[[0, 0]] = 1e-5;
        bm[[1, 0]] = 1e-5;
        bm[[2, 1]] = 1e-5;

        let mut _c_e: Array2<f32>;

        let pb = ProgressBar::new(m1_n_mode as u64);
        pb.set_style(ProgressStyle::default_bar().template("{msg}"));
        for k_step in 0..20 {
            let msg = String::from("Step #") + &k_step.to_string();
            pb.println(msg);

            for sid in 1..8 {
                //print!("{}", sid);••••••••••••
                id = sid - 1;
                t_xyz[0] = rbm[[id, 0]] as f64;
                t_xyz[1] = rbm[[id, 1]] as f64;
                t_xyz[2] = rbm[[id, 2]] as f64;
                r_xyz[0] = rbm[[id, 3]] as f64;
                r_xyz[1] = rbm[[id, 4]] as f64;
                r_xyz[2] = rbm[[id, 5]] as f64;
                gmt.set_m1_segment_state(sid as i32, &t_xyz, &r_xyz);
                if m1_n_mode > 0 {
                    for k_bm in 0..m1_n_mode {
                        let idx = id * m1_n_mode as usize + k_bm as usize;
                        a[idx as usize] = bm[[id, k_bm as usize]] as f64;
                    }
                }
                id += 7;
                t_xyz[0] = rbm[[id, 0]] as f64;
                t_xyz[1] = rbm[[id, 1]] as f64;
                t_xyz[2] = rbm[[id, 2]] as f64;
                r_xyz[0] = rbm[[id, 3]] as f64;
                r_xyz[1] = rbm[[id, 4]] as f64;
                r_xyz[2] = rbm[[id, 5]] as f64;
                gmt.set_m2_segment_state(sid as i32, &t_xyz, &r_xyz);
            }
            gmt.set_m1_modes(&mut a);

            src_wfe_rms = src.through(&mut gmt).wfe_rms_10e(-9)[0];
            onaxis_pssn = src_pssn.reset(&mut src);

            //println!(
            //   "WFE RMS: {:.3} nm; PSSn: {:.5}",
            //    src_wfe_rms * 1e9,
            //    src_pssn.eval()
            //);
            let msg = String::from("WFE RMS: ")
                + &src_wfe_rms.to_string()
                + "nm; PSSn: "
                + &onaxis_pssn.to_string();
            pb.println(msg);

            gs.through(&mut gmt).through(&mut gwfs);
            /*
            for _k_atm in (0..atm_sample).progress() {
                gs.wavefront.reset();
                gs.opd2phase();
                atm.rayTracing1(&mut gs, px, pupil_sampling, px, pupil_sampling, atm_time);
                wfs.propagate(&mut gs);
                atm_time += atm_sampling;
            }*/
            gwfs.process();
            let slopes = Array::from_shape_vec((n_c, 1), gwfs.centroids.clone()).unwrap();

            _c_e = __m.dot(&slopes);
            let q = gain
                * _c_e
                    .slice(s![..84, ..])
                    .to_owned()
                    .into_shape((14, 6))
                    .unwrap();
            rbm -= &q.view();
            if m1_n_mode > 0 {
                let q = gain
                    * _c_e
                        .slice(s![84.., ..])
                        .to_owned()
                        .into_shape((7, m1_n_mode as usize))
                        .unwrap();
                bm -= &q.view();
            }
        }
    }
}
