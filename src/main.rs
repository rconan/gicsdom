#[macro_use]
extern crate crossbeam_channel;

use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::ceo;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{s, stack, Array, Array2, Axis};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::f32;
use std::thread;
use std::time::{Duration, Instant};

#[derive(Debug)]
struct GmtState {
    rbm: Array2<f32>,
    bm: Array2<f32>,
}

fn main() {
    // Create a zero-capacity channel.
    let plant_to_sh48: (Sender<GmtState>, Receiver<GmtState>) = bounded(0);
    let sh48_plan_chat: (Sender<GmtState>, Receiver<GmtState>) = bounded(0);
    let plant_to_science: (Sender<GmtState>, Receiver<GmtState>) = unbounded();
    let science_plan_chat = bounded(0);

    // SH48
    let r = plant_to_sh48.1;
    let s = sh48_plan_chat.0;
    thread::spawn(move || {
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

        // ---------------------------------------------------
        // SYSTEM CALIBRATION
        let sh48_from_calibrate = bounded(0);
        let sp = sh48_from_calibrate.0;
        thread::spawn(move || {
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

            //println!("WFS centroids #     : {}", gwfs.n_centroids);
            //println!("WFS centroids length: {}", gwfs.centroids.len());
            //println!("CCM length: {}", c_c_m.len());

            let _d_m = Array2::from_shape_vec((c_c_m.len() / n_c, n_c), c_c_m)
                .unwrap()
                .t()
                .to_owned()
                * 1e6;
            let mut _d = (_d_p - _d_m) * 0.5;
            println!(" in {}s", now.elapsed().as_secs());

            //println!("{:?}", _d);
            println!("shape={:?}, strides={:?}", _d.shape(), _d.strides());
            //        println!("d sum: {}",_d.into.sum());

            print!("SVD decomposition");
            let now = Instant::now();
            let n_sv = n_rbm + m1_n_mode as usize * 7;
            let (u, sig, v_t) = _d.svddc_inplace(UVTFlag::Some).unwrap();
            println!(" in {}ms", now.elapsed().as_millis());
            //println!("Singular values:\n{}", sig);
            let mut i_sig = sig.mapv(|x| 1.0 / x);
            for k in 0..14 {
                i_sig[n_sv - k - 1] = 0.0;
            }

            let _u = u.unwrap();
            let _vt = v_t.unwrap();

            let l_sv = Array2::from_diag(&i_sig);
            let z: Array2<f32> = Array2::zeros((n_sv, n_c - n_sv));
            let ss = stack(Axis(1), &[l_sv.view(), z.view()]).unwrap();
            //println!("SS shape: {:?}", ss.shape());
            print!("Computing the pseudo-inverse");
            let now = Instant::now();
            let __m: Array2<f32> = _vt.t().dot(&l_sv.dot(&_u.t()));
            println!(" in {}ms", now.elapsed().as_millis());

            sp.send(__m).unwrap();
        });
        // ---------------------------------------------------

        let rp = sh48_from_calibrate.1;
        //println!("{:?}", rp.recv());
        let __m = rp.recv().unwrap();

        //println!("{:?}", r.recv());

        let n_c: usize = (n_side_lenslet.pow(2) as usize) * 2 * n_gs as usize;
        let gain = 0.5f32;
        let mut _c_e: Array2<f32>;
        let mut t_xyz = vec![0.0; 3];
        let mut r_xyz = vec![0.0; 3];
        let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
        let mut id = 0;
        let mut rbm = Array2::<f32>::zeros((14, 6));
        let mut bm = Array2::<f32>::zeros((7, m1_n_mode as usize));

        loop {

            let gstate = match r.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            for sid in 1..8 {
                //print!("{}", sid);••••••••••••
                id = sid - 1;
                t_xyz[0] = gstate.rbm[[id, 0]] as f64;
                t_xyz[1] = gstate.rbm[[id, 1]] as f64;
                t_xyz[2] = gstate.rbm[[id, 2]] as f64;
                r_xyz[0] = gstate.rbm[[id, 3]] as f64;
                r_xyz[1] = gstate.rbm[[id, 4]] as f64;
                r_xyz[2] = gstate.rbm[[id, 5]] as f64;
                gmt.set_m1_segment_state(sid as i32, &t_xyz, &r_xyz);
                if m1_n_mode > 0 {
                    for k_bm in 0..m1_n_mode {
                        let idx = id * m1_n_mode as usize + k_bm as usize;
                        a[idx as usize] = gstate.bm[[id, k_bm as usize]] as f64;
                    }
                }
                id += 7;
                t_xyz[0] = gstate.rbm[[id, 0]] as f64;
                t_xyz[1] = gstate.rbm[[id, 1]] as f64;
                t_xyz[2] = gstate.rbm[[id, 2]] as f64;
                r_xyz[0] = gstate.rbm[[id, 3]] as f64;
                r_xyz[1] = gstate.rbm[[id, 4]] as f64;
                r_xyz[2] = gstate.rbm[[id, 5]] as f64;
                gmt.set_m2_segment_state(sid as i32, &t_xyz, &r_xyz);
            }
            gmt.set_m1_modes(&mut a);

            gs.through(&mut gmt).through(&mut gwfs);
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

            let mut gstate = GmtState {
                rbm: Array2::<f32>::zeros((14, 6)),
                bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
            };
            gstate.rbm = rbm.clone();
            gstate.bm = bm.clone();
            s.send(gstate).unwrap();
 
        }

        drop(gmt);
        drop(gwfs);
        drop(gs);

        let mut gstate = GmtState {
            rbm: Array2::<f32>::zeros((14, 6)),
            bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
        };
        s.send(gstate).unwrap();
    });

    // SCIENCE
    let r = plant_to_science.1;
    let s = science_plan_chat.0;
    thread::spawn(move || {
        let pupil_size = 25.5;
        let pupil_sampling = 512;
        let m1_n_mode = 0;

        // GMT initialization
        let mut gmt = ceo::Gmt::new(m1_n_mode, None);
        gmt.build();

        // Science initialization
        let mut src = ceo::Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0.0], vec![0.0], vec![0.0]);

        /*
        for k in 0..10 {
            thread::sleep(Duration::from_secs(5));
            println!("WFE RMS: {}", src.through(&mut gmt).wfe_rms_10e(-9)[0]);
        }
        */

        let mut src_wfe_rms = 0f32;
        src_wfe_rms = src.through(&mut gmt).wfe_rms_10e(-9)[0];
        let mut onaxis_pssn = 0f32;
        let mut src_pssn = ceo::PSSn::new(15e-2, 25.0, 0.0);
        src_pssn.build(&mut src);
        println!(
            "WFE RMS: {:.3} nm; PSSn: {}",
            src_wfe_rms,
            src_pssn.reset(&mut src)
        );

        //        println!("{:?}", r.recv());
        loop {
            let mut t_xyz = vec![0.0; 3];
            let mut r_xyz = vec![0.0; 3];
            let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
            let mut id = 0;

            let gstate = match r.recv() {
                Ok(state) => state,
                Err(error) => break,
            };

            for sid in 1..8 {
                //print!("{}", sid);••••••••••••
                id = sid - 1;
                t_xyz[0] = gstate.rbm[[id, 0]] as f64;
                t_xyz[1] = gstate.rbm[[id, 1]] as f64;
                t_xyz[2] = gstate.rbm[[id, 2]] as f64;
                r_xyz[0] = gstate.rbm[[id, 3]] as f64;
                r_xyz[1] = gstate.rbm[[id, 4]] as f64;
                r_xyz[2] = gstate.rbm[[id, 5]] as f64;
                gmt.set_m1_segment_state(sid as i32, &t_xyz, &r_xyz);
                if m1_n_mode > 0 {
                    for k_bm in 0..m1_n_mode {
                        let idx = id * m1_n_mode as usize + k_bm as usize;
                        a[idx as usize] = gstate.bm[[id, k_bm as usize]] as f64;
                    }
                }
                id += 7;
                t_xyz[0] = gstate.rbm[[id, 0]] as f64;
                t_xyz[1] = gstate.rbm[[id, 1]] as f64;
                t_xyz[2] = gstate.rbm[[id, 2]] as f64;
                r_xyz[0] = gstate.rbm[[id, 3]] as f64;
                r_xyz[1] = gstate.rbm[[id, 4]] as f64;
                r_xyz[2] = gstate.rbm[[id, 5]] as f64;
                gmt.set_m2_segment_state(sid as i32, &t_xyz, &r_xyz);
            }
            gmt.set_m1_modes(&mut a);

            src_wfe_rms = src.through(&mut gmt).wfe_rms_10e(-9)[0];
            onaxis_pssn = src_pssn.reset(&mut src);

            println!(
                "@science => WFE RMS: {}nm; PSSn: {}",
                &src_wfe_rms.to_string(),
                &onaxis_pssn.to_string()
            );
            thread::sleep(Duration::from_secs(1));
        }
        //        s.send("Hi from SCIENCE!").unwrap();

        drop(gmt);
        drop(src);

        s.send("Bye from SCIENCE!".to_string()).unwrap();
    });

    // PLANT

    // GMT initial conditions
    let m1_n_mode = 27;
    let mut gstate0 = GmtState {
        rbm: Array2::<f32>::zeros((14, 6)),
        bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
    };
    gstate0.rbm[(0, 0)] = 1e-5;
    gstate0.rbm[(8, 4)] = 1e-6;
    gstate0.rbm[(2, 3)] = -1e-6;
    gstate0.rbm[(11, 1)] = -1e-5;
    gstate0.bm[[0, 0]] = 1e-5;
    gstate0.bm[[1, 0]] = 1e-5;
    gstate0.bm[[2, 1]] = 1e-5;
    
    let s_science = plant_to_science.0;
    let s_sh48 = plant_to_sh48.0;
    let mut gstate = GmtState {
        rbm: Array2::<f32>::zeros((14, 6)),
        bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
    };
    gstate.rbm = gstate0.rbm.clone();
    gstate.bm = gstate0.bm.clone();

    let r = sh48_plan_chat.1;
     for _ in 0..20 {

        let mut gstate_science = GmtState {
            rbm: Array2::<f32>::zeros((14, 6)),
            bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
        };
        gstate_science.rbm = gstate.rbm.clone();
        gstate_science.bm = gstate.bm.clone();
        s_science.send(gstate_science).unwrap();

        let mut gstate_sh48 = GmtState {
            rbm: Array2::<f32>::zeros((14, 6)),
            bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
        };
        gstate_sh48.rbm = gstate.rbm.clone();
        gstate_sh48.bm = gstate.bm.clone();
        s_sh48.send(gstate_sh48).unwrap();

        let new_gstate = match r.recv() {
            Ok(state) => state,
            Err(_) => break,
        };
        gstate.rbm = gstate0.rbm.clone() + new_gstate.rbm;
        gstate.bm = gstate0.bm.clone() + new_gstate.bm;

    }

    drop(s_science);
    drop(s_sh48);

    r.recv().unwrap();
    let r = science_plan_chat.1;
    println!("{:?}", r.recv());
}
