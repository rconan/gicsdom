//extern crate ndarray;
//extern crate ndarray_linalg;

use std::ffi::CString;
use std::time::Instant;
use std::{f32, mem};
//use gicsdom::optics_path_gsh48::OpticsPathGSH48;
//use gicsdom::optics_path_sh48::OpticsPathSH48;
use gicsdom::{
    atmosphere, dev2host, geometricShackHartmann, gmt_m1, gmt_m2, modes, pssn, shackHartmann,
    source, vector, zernikeS,
};
//use ndarray::prelude::*;
use ndarray::{s, stack, Array, Array2, Axis};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
//use rand::distributions::{Normal, Uniform};
//use rand::thread_rng;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};

fn main() {
    unsafe {
        //        let rng = thread_rng();

        let n_side_lenslet = 48;
        let n_px_lenslet = 16;
        //let n_px = n_side_lenslet*16 + 1;
        let pupil_size = 25.5;
        let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;

        // M1 initialization
        println!("Initializing M1!");
        let m1_n_mode: usize = 27;
        let mode_type = CString::new("bending modes").unwrap();
        let mut m1_modes: modes = mem::zeroed();
        let mut m1: gmt_m1 = mem::zeroed();
        //        let mut a: Vec<f64> = vec![0.0];
        m1_modes //.setup1(0, a.as_mut_ptr(), 7);
            .setup(mode_type.into_raw() as *mut i8, 7, m1_n_mode as i32);
        m1.setup1(&mut m1_modes);

        // M2 initialization
        println!("Initializing M2!");
        let mut m2_modes: zernikeS = mem::zeroed();
        let mut m2: gmt_m2 = mem::zeroed();
        let mut a: Vec<f64> = vec![0.0];
        m2_modes.setup1(0, a.as_mut_ptr(), 7);
        m2.setup2(&mut m2_modes);

        // Guide star initialization
        println!("Initializing guide star!");
        let origin = vector {
            x: 0.0,
            y: 0.0,
            z: 25.0,
        };
        let band = CString::new("V").unwrap();
        let mut magnitude: Vec<f32> = vec![0.0, 0.0, 0.0];
        let z = 6. * f32::consts::PI / 180. / 60.;
        let mut zenith: Vec<f32> = vec![z, z, z];
        let a = 2. * f32::consts::PI / 3.;
        let mut azimuth: Vec<f32> = vec![0.0, a, 2. * a];
        let n_gs = zenith.len() as i32;
        let mut gs: source = mem::zeroed();
        gs.setup7(
            band.into_raw(),
            magnitude.as_mut_ptr(),
            zenith.as_mut_ptr(),
            azimuth.as_mut_ptr(),
            f32::INFINITY,
            n_gs,
            pupil_size,
            pupil_sampling,
            origin,
        );

        // GWFS initialization
        println!("Initializing wavefront sensor!");
        let wfs_exposure_time = 30f32;
        let d: f32 = pupil_size as f32 / n_side_lenslet as f32;
        println!("d={}", d);
        let mut gwfs: geometricShackHartmann = mem::zeroed();
        gwfs.setup(n_side_lenslet, d, n_gs);
        println!("Calibration wavefront sensor!");
        // GWFS calibration
        gs.wavefront.reset();
        gs.reset_rays();
        m2.blocking(&mut gs.rays);
        m1.trace(&mut gs.rays);
        //gs.rays.gmt_truss_onaxis();
        m2.trace(&mut gs.rays);
        gs.rays.to_sphere1(-5.830, 2.197173);
        gs.opd2phase();
        gwfs.calibrate(&mut gs, 0.0);
        gwfs.reset();

        // Atmosphere initialization
        println!("Initializing atmosphere!");
        let r_not = 15e-2;
        //let r0_time_series[15];
        let l_not = 30f32;
        let screen_size = 26f32;
        let screen_sampling = 346;
        let field_size = 20.0 * f32::consts::PI / 180.0 / 60.0;
        let duration = 15f32;
        let n_duration = 4;
        let pathtoscreens = CString::new("/home/rconan/DATA/gmtAtmosphereL030_1.bin").unwrap();
        //char pathtoscreens[] = "/home/ubuntu/DATA/gmtAtmosphereL030_1579812935.bin";
        let mut atm_time = 0.0;
        let atm_sampling = 1e-2;
        let atm_sample = (wfs_exposure_time / atm_sampling) as u64;
        let mut atm: atmosphere = mem::zeroed();
        atm.gmt_setup3(
            r_not,
            l_not,
            screen_size,
            screen_sampling,
            field_size,
            duration,
            pathtoscreens.into_raw(),
            n_duration,
        );

        // ----------------------------------------------------------------------------
        // GMT CALIBRATION
        println!("Calibrating the GMT ...");
        let now = Instant::now();
        let n_rbm: usize = 84;
        let mut c_c_p: Vec<f32> = Vec::new();
        let mut c_c_m: Vec<f32> = Vec::new();
        let n_c: usize = (n_side_lenslet.pow(2) as usize) * 2 * n_gs as usize;
        m1.reset();
        m2.reset();
        let pb = ProgressBar::new(7);
        pb.set_style(ProgressStyle::default_bar().template("{msg}"));
        println!("Rigid body motion ...");
        if n_rbm > 0 {
            for mid in 1..3 {
                let m_ = String::from("M") + &mid.to_string();
                pb.set_message(&m_);
                for sid in (1..8).progress() {
                    for tr in 1..3 {
                        for a in 0..3 {
                            m1.reset();
                            m2.reset();
                            let mut t_xyz = vector {
                                x: 0.0,
                                y: 0.0,
                                z: 0.0,
                            };
                            let mut r_xyz = vector {
                                x: 0.0,
                                y: 0.0,
                                z: 0.0,
                            };
                            if tr == 1 {
                                if a == 0 {
                                    t_xyz.x = 1e-6;
                                }
                                if a == 1 {
                                    t_xyz.y = 1e-6;
                                }
                                if a == 2 {
                                    t_xyz.z = 1e-6;
                                }
                            }
                            if tr == 2 {
                                if a == 0 {
                                    r_xyz.x = 1e-6;
                                }
                                if a == 1 {
                                    r_xyz.y = 1e-6;
                                }
                                if a == 2 {
                                    r_xyz.z = 1e-6;
                                }
                            }
                            if mid == 1 {
                                m1.update(t_xyz, r_xyz, sid);
                            }
                            if mid == 2 {
                                m2.update(t_xyz, r_xyz, sid);
                            }
                            gs.wavefront.reset();
                            gs.reset_rays();
                            m2.blocking(&mut gs.rays);
                            m1.trace(&mut gs.rays);
                            // gs.rays.gmt_truss_onaxis();
                            m2.trace(&mut gs.rays);
                            gs.rays.to_sphere1(-5.830, 2.197173);
                            gs.opd2phase();
                            gwfs.propagate(&mut gs);
                            gwfs.process();
                            let mut c: Vec<f32> = vec![0.0; n_c];
                            dev2host(c.as_mut_ptr(), gwfs.data_proc.d__c, n_c as i32);
                            gwfs.reset();
                            //let sum1: f32 = c.iter().sum();
                            //println!("Centroids sum: {:.16}", sum1);
                            c_c_p.append(&mut c);
                        }
                    }
                }
            }
        }

        if m1_n_mode > 0 {
            println!("Bending modes ...");
            let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
            m1.reset();
            m2.reset();
            let pb = ProgressBar::new(m1_n_mode as u64);
            pb.set_style(ProgressStyle::default_bar().template("{msg}"));
            for sid in 0..7 {
                let s_ = String::from("S") + &(sid + 1).to_string();
                pb.set_message(&s_);
                for k_a in (0..m1_n_mode).progress() {
                    let k = k_a + sid * m1_n_mode;
                    a[k] = 1e-6;
                    m1_modes.update(a.as_mut_ptr());
                    gs.wavefront.reset();
                    gs.reset_rays();
                    m2.blocking(&mut gs.rays);
                    m1.trace(&mut gs.rays);
                    // gs.rays.gmt_truss_onaxis();
                    m2.trace(&mut gs.rays);
                    gs.rays.to_sphere1(-5.830, 2.197173);
                    gs.opd2phase();
                    gwfs.propagate(&mut gs);
                    gwfs.process();
                    let mut c: Vec<f32> = vec![0.0; n_c];
                    dev2host(c.as_mut_ptr(), gwfs.data_proc.d__c, n_c as i32);
                    gwfs.reset();
                    a[k] = 0.0;
                    m1_modes.update(a.as_mut_ptr());
                    &c_c_p.append(&mut c);
                }
            }
        }

        let _d_p = Array2::from_shape_vec((c_c_p.len() / n_c, n_c), c_c_p)
            .unwrap()
            .t()
            .to_owned()
            * 1e6;

        println!("Rigid body motion ...");
        let pb = ProgressBar::new(7);
        pb.set_style(ProgressStyle::default_bar().template("{msg}"));
        if n_rbm > 0 {
            for mid in 1..3 {
                let m_ = String::from("M") + &mid.to_string();
                pb.set_message(&m_);
                for sid in (1..8).progress() {
                    for tr in 1..3 {
                        for a in 0..3 {
                            m1.reset();
                            m2.reset();
                            let mut t_xyz = vector {
                                x: 0.0,
                                y: 0.0,
                                z: 0.0,
                            };
                            let mut r_xyz = vector {
                                x: 0.0,
                                y: 0.0,
                                z: 0.0,
                            };
                            if tr == 1 {
                                if a == 0 {
                                    t_xyz.x = -1e-6;
                                }
                                if a == 1 {
                                    t_xyz.y = -1e-6;
                                }
                                if a == 2 {
                                    t_xyz.z = -1e-6;
                                }
                            }
                            if tr == 2 {
                                if a == 0 {
                                    r_xyz.x = -1e-6;
                                }
                                if a == 1 {
                                    r_xyz.y = -1e-6;
                                }
                                if a == 2 {
                                    r_xyz.z = -1e-6;
                                }
                            }
                            if mid == 1 {
                                m1.update(t_xyz, r_xyz, sid);
                            }
                            if mid == 2 {
                                m2.update(t_xyz, r_xyz, sid);
                            }
                            gs.wavefront.reset();
                            gs.reset_rays();
                            m2.blocking(&mut gs.rays);
                            m1.trace(&mut gs.rays);
                            // gs.rays.gmt_truss_onaxis();
                            m2.trace(&mut gs.rays);
                            gs.rays.to_sphere1(-5.830, 2.197173);
                            gs.opd2phase();
                            gwfs.propagate(&mut gs);
                            gwfs.process();
                            let mut c: Vec<f32> = vec![0.0; n_c];
                            dev2host(c.as_mut_ptr(), gwfs.data_proc.d__c, n_c as i32);
                            gwfs.reset();
                            //let sum2: f32 = c.iter().sum();
                            //println!("Centroids sum: {:.16}", sum2);
                            //println!("Diff. centroids sum: {:.14}", sum2 - sum1);
                            c_c_m.append(&mut c);
                        }
                    }
                }
            }
        }

        if m1_n_mode > 0 {
            println!("Bending modes ...");
            let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
            m1.reset();
            m2.reset();
            let pb = ProgressBar::new(m1_n_mode as u64);
            pb.set_style(ProgressStyle::default_bar().template("{msg}"));
            for sid in 0..7 {
                let s_ = String::from("S") + &(sid + 1).to_string();
                pb.set_message(&s_);
                for k_a in (0..m1_n_mode).progress() {
                    let k = k_a + sid * m1_n_mode;
                    a[k] = -1e-6;
                    m1_modes.update(a.as_mut_ptr());
                    gs.wavefront.reset();
                    gs.reset_rays();
                    m2.blocking(&mut gs.rays);
                    m1.trace(&mut gs.rays);
                    // gs.rays.gmt_truss_onaxis();
                    m2.trace(&mut gs.rays);
                    gs.rays.to_sphere1(-5.830, 2.197173);
                    gs.opd2phase();
                    gwfs.propagate(&mut gs);
                    gwfs.process();
                    let mut c: Vec<f32> = vec![0.0; n_c];
                    dev2host(c.as_mut_ptr(), gwfs.data_proc.d__c, n_c as i32);
                    gwfs.reset();
                    a[k] = 0.0;
                    m1_modes.update(a.as_mut_ptr());
                    &c_c_m.append(&mut c);
                }
            }
        }

        let _d_m = Array2::from_shape_vec((c_c_m.len() / n_c, n_c), c_c_m)
            .unwrap()
            .t()
            .to_owned()
            * 1e6;
        let mut _d = (_d_p - _d_m) * 0.5;
        println!(" in {}s", now.elapsed().as_secs());

        println!("{:?}", _d);
        println!("shape={:?}, strides={:?}", _d.shape(), _d.strides());

        print!("SVD decomposition");
        let now = Instant::now();
        let n_sv = n_rbm + m1_n_mode * 7;
        let (u, sig, v_t) = _d.svddc_inplace(UVTFlag::Some).unwrap();
        println!(" in {}ms", now.elapsed().as_millis());
        //        println!("Singular values:\n{}", sig);
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
        let mut t_xyz = vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut r_xyz = vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];
        m1.reset();
        m2.reset();
        t_xyz.y = 1e-6;
        m1.update(t_xyz, r_xyz, 1);
        r_xyz.x = 1e-6;
        m2.update(t_xyz, r_xyz, 2);
        if m1_n_mode > 0 {
            a[0] = 1e-6;
            a[5] = 1e-6;
            a[11] = 1e-6;
            m1_modes.update(a.as_mut_ptr());
        }
        gs.wavefront.reset();
        gs.reset_rays();
        m2.blocking(&mut gs.rays);
        m1.trace(&mut gs.rays);
        // gs.rays.gmt_truss_onaxis();
        m2.trace(&mut gs.rays);
        gs.rays.to_sphere1(-5.830, 2.197173);
        gs.opd2phase();
        gwfs.propagate(&mut gs);
        gwfs.process();
        let mut c: Vec<f32> = vec![0.0; n_c];
        dev2host(c.as_mut_ptr(), gwfs.data_proc.d__c, n_c as i32);
        gwfs.reset();

        let slopes = Array::from_shape_vec((n_c, 1), c).unwrap();
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
                bm.to_owned().into_shape((7, m1_n_mode)).unwrap() * 1e6
            );
        }
        // ------------------------------------------------------------------------------------------------------
        let closed_loop = true;
        if closed_loop {
            //        let normal: Normal = Normal::new(0,1e-7);
            //        let mut com = MatrixMN::<f64, U6, U14>::from_distribution(&normal,&mut rng);
            //let mut com: Array2<f32> = Array2::zeros((6,14));
            let gain = 0.5f32;

            // Science initialization
            let band = CString::new("V").unwrap();
            let mut magnitude: Vec<f32> = vec![0.0];
            let mut zenith: Vec<f32> = vec![0.0];
            let mut azimuth: Vec<f32> = vec![0.0];
            let n_gs = zenith.len() as i32;
            let mut src: source = mem::zeroed();
            src.setup7(
                band.into_raw(),
                magnitude.as_mut_ptr(),
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                f32::INFINITY,
                n_gs,
                pupil_size,
                pupil_sampling,
                origin,
            );
            let mut src_wfe_rms = 0f32;
            let mut src_pssn: pssn = mem::zeroed();

            let mut t = vector {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            };
            let mut r = vector {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            };
            let mut a: Vec<f64> = vec![0.0; 7 * m1_n_mode as usize];

            m1.reset();
            m2.reset();
            m1_modes.update(a.as_mut_ptr());
            src.wavefront.reset();
            src.reset_rays();
            m2.blocking(&mut src.rays);
            m1.trace(&mut src.rays);
            //src.rays.gmt_truss_onaxis();
            m2.trace(&mut src.rays);
            src.rays.to_sphere1(-5.830, 2.197173);
            src.opd2phase();
            src.wavefront.rms(&mut src_wfe_rms);

            src_pssn.setup(&mut src, 15e-2, 25.0);
            src_pssn.N_O = 0;
            src_pssn.otf(&mut src);
            println!(
                "WFE RMS: {:.3} nm; PSSn: {}",
                src_wfe_rms * 1e9,
                src_pssn.eval()
            );

            // WFS initialization
            println!("Initializing wavefront sensor!");
            let mut wfs: shackHartmann = mem::zeroed();
            let d: f32 = pupil_size as f32 / n_side_lenslet as f32;
            wfs.setup(n_side_lenslet, n_px_lenslet, d, 2, 24, 3, 3);
            gs.wavefront.reset();
            println!("Calibration wavefront sensor!");
            // WFS calibration
            gs.reset_rays();
            m2.blocking(&mut gs.rays);
            m1.trace(&mut gs.rays);
            //gs.rays.gmt_truss_onaxis();
            m2.trace(&mut gs.rays);
            gs.rays.to_sphere1(-5.830, 2.197173);
            gs.opd2phase();
            gs.fwhm = 3.16;
            wfs.calibrate(&mut gs, 0.0);
            wfs.camera.reset();

            //println!("{:.2}", com[(0,0)] * 1e6);
            let mut id = 0;
            let px = (pupil_size/(pupil_sampling as f64 - 1.0)) as f32;

            // GMT initial conditions
            let mut rbm = Array2::<f32>::zeros((14, 6));
            rbm[(0, 0)] = 1e-5;
            rbm[(8, 4)] = 1e-6;
            rbm[(2, 3)] = -1e-6;
            rbm[(11, 1)] = -1e-5;
            let mut bm = Array2::<f32>::zeros((7, m1_n_mode));
            bm[[0, 0]] = 1e-5;
            bm[[1, 0]] = 1e-5;
            bm[[2, 1]] = 1e-5;

            let mut slopes = Array2::<f32>::zeros((n_c, 1));
            let mut _c_e: Array2<f32>;

            let pb = ProgressBar::new(m1_n_mode as u64);
            pb.set_style(ProgressStyle::default_bar().template("{msg}"));
            for k_step in 0..2 {
                let msg = String::from("Step #") + &k_step.to_string();
                pb.println(msg);

                for sid in 1..8 {
                    //print!("{}", sid);
                    id = sid - 1;
                    t.x = rbm[[id, 0]] as f64;
                    t.y = rbm[[id, 1]] as f64;
                    t.z = rbm[[id, 2]] as f64;
                    r.x = rbm[[id, 3]] as f64;
                    r.y = rbm[[id, 4]] as f64;
                    r.z = rbm[[id, 5]] as f64;
                    m1.update(t, r, sid as i32);
                    if m1_n_mode > 0 {
                        for k_bm in 0..m1_n_mode {
                            a[id * m1_n_mode + k_bm] = bm[[id, k_bm]] as f64;
                        }
                    }
                    id += 7;
                    t.x = rbm[[id, 0]] as f64;
                    t.y = rbm[[id, 1]] as f64;
                    t.z = rbm[[id, 2]] as f64;
                    r.x = rbm[[id, 3]] as f64;
                    r.y = rbm[[id, 4]] as f64;
                    r.z = rbm[[id, 5]] as f64;
                    m2.update(t, r, sid as i32);
                }
                m1_modes.update(a.as_mut_ptr());
                gs.wavefront.reset();

                src.wavefront.reset();
                src.reset_rays();
                m2.blocking(&mut src.rays);
                m1.trace(&mut src.rays);
                //src.rays.gmt_truss_onaxis();
                m2.trace(&mut src.rays);
                src.rays.to_sphere1(-5.830, 2.197173);
                src.opd2phase();
                src.wavefront.rms(&mut src_wfe_rms);
                src_pssn.N_O = 0;
                src_pssn.otf(&mut src);
                //println!(
                //   "WFE RMS: {:.3} nm; PSSn: {:.5}",
                //    src_wfe_rms * 1e9,
                //    src_pssn.eval()
                //);
                let msg = String::from("WFE RMS: ") + &(src_wfe_rms * 1e9).to_string() + "nm; PSSn: " + &src_pssn.eval().to_string();
                pb.println(msg);

                gs.reset_rays();
                m2.blocking(&mut gs.rays);
                m1.trace(&mut gs.rays);
                // gs.rays.gmt_truss_onaxis();
                m2.trace(&mut gs.rays);
                gs.rays.to_sphere1(-5.830, 2.197173);
                for _k_atm in (0..atm_sample).progress() {
                    gs.wavefront.reset();
                    gs.opd2phase();
                    atm.rayTracing1(&mut gs,px,pupil_sampling,px,pupil_sampling,atm_time);
                    wfs.propagate(&mut gs);
                    atm_time += atm_sampling;
                }
                wfs.process();
                dev2host(slopes.as_mut_ptr(), wfs.data_proc.d__c, n_c as i32);
                wfs.camera.reset();

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
                            .into_shape((7, m1_n_mode))
                            .unwrap();
                    bm -= &q.view();
                }
            }
        }

        m1_modes.cleanup();
        m1.cleanup();
        m2_modes.cleanup();
        m2.cleanup();
        gs.cleanup();
        gwfs.cleanup();
        atm.cleanup();
    }
}
