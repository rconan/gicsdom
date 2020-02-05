extern crate nalgebra as na;
extern crate ndarray;
extern crate ndarray_linalg;

use std::ffi::CString;
use std::{f32, mem};
use std::time::{Duration, Instant};
//use gicsdom::optics_path_gsh48::OpticsPathGSH48;
//use gicsdom::optics_path_sh48::OpticsPathSH48;
use gicsdom::{
    dev2host, geometricShackHartmann, gmt_m1, gmt_m2, modes, pssn, source, vector, zernikeS,
};
//use ndarray::prelude::*;
use ndarray::{stack,Axis, Array, Array2};
use ndarray_linalg::*;
//use rand::distributions::{Normal, Uniform};
//use rand::thread_rng;

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
        let m1_n_mode = 27;
        let mode_type = CString::new("bending modes").unwrap();
        let mut m1_modes: zernikeS = mem::zeroed();
        let mut m1: gmt_m1 = mem::zeroed();
        let mut a: Vec<f64> = vec![0.0];
        m1_modes.setup1(0, a.as_mut_ptr(), 7);
        //    .setup(mode_type.into_raw() as *mut i8, 7, m1_n_mode);
        m1.setup2(&mut m1_modes);

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

        // WFS initialization
        println!("Initializing wavefront sensor!");
        let d: f32 = pupil_size as f32 / n_side_lenslet as f32;
        println!("d={}", d);
        let mut wfs: geometricShackHartmann = mem::zeroed();
        wfs.setup(n_side_lenslet, d, n_gs);
        println!("Calibration wavefront sensor!");
        // WFS calibration
        gs.wavefront.reset();
        gs.reset_rays();
        m2.blocking(&mut gs.rays);
        m1.trace(&mut gs.rays);
        //gs.rays.gmt_truss_onaxis();
        m2.trace(&mut gs.rays);
        gs.rays.to_sphere1(-5.830, 2.197173);
        gs.opd2phase();
        wfs.calibrate(&mut gs, 0.0);
        wfs.reset();

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
        let mut c_c_p: Vec<f32> = Vec::new();
        let mut c_c_m: Vec<f32> = Vec::new();

        let n_c: usize = (n_side_lenslet.pow(2) as usize) * 2 * n_gs as usize;
        m1.reset();
        m2.reset();
        for mid in 1..3 {
            for sid in 1..8 {
                print!("{}", sid);
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
                        wfs.propagate(&mut gs);
                        wfs.process();
                        let mut c: Vec<f32> = vec![0.0; n_c];
                        dev2host(c.as_mut_ptr(), wfs.data_proc.d__c, n_c as i32);
                        wfs.reset();
                        //let sum1: f32 = c.iter().sum();
                        //println!("Centroids sum: {:.16}", sum1);
                        c_c_p.append(&mut c);
                    }
                }
            }
        }
        println!("");
        let _d_p = Array2::from_shape_vec((c_c_p.len() / n_c, n_c), c_c_p)
            .unwrap()
            .t()
            .to_owned()
            * 1e6;

        for mid in 1..3 {
            for sid in 1..8 {
                print!("{}", sid);
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
                        wfs.propagate(&mut gs);
                        wfs.process();
                        let mut c: Vec<f32> = vec![0.0; n_c];
                        dev2host(c.as_mut_ptr(), wfs.data_proc.d__c, n_c as i32);
                        wfs.reset();
                        //let sum2: f32 = c.iter().sum();
                        //println!("Centroids sum: {:.16}", sum2);
                        //println!("Diff. centroids sum: {:.14}", sum2 - sum1);
                        c_c_m.append(&mut c);
                    }
                }
            }
        }
        println!("");
        let _d_m = Array::from_shape_vec((c_c_m.len() / n_c, n_c), c_c_m)
            .unwrap()
            .t()
            .to_owned()
            * 1e6;

        let _d = (_d_p - _d_m) * 0.5;
        //println!("D shape: {:?}",_d.shape());
        print!("SVD decomposition");
        let now = Instant::now();
        let (u, sig, v_t) = _d.svd(true, true).unwrap();
        println!("in {}s", now.elapsed().as_secs());
        //println!("U shape: {:?}",u.unwrap().shape());
        //println!("V.T shape: {:?}",v_t.unwrap().shape());
        //println!("S shape: {:?}",sig.shape());
        println!("Singular values:\n{}", sig);
        //let _u = u.unwrap();
        //let _v_t = v_t.unwrap();
        //let _m = _u.dot(&sig.dot(&_v_t));

        let mut i_sig = sig.mapv(|x| 1.0 / x);
        for k in 0..14 {
            i_sig[83 - k] = 0.0;
        }
        //println!("Singular values:\n{}", i_sig);
        let s = Array2::from_diag(&i_sig);
        let z: Array2<f32> = Array2::zeros((84, n_c - 84));
        let ss = stack(Axis(1), &[s.view(), z.view()]).unwrap();
        println!("SS shape: {:?}", ss.shape());
        print!("Computing the pseudo-inverse ");
        let now = Instant::now();
        let __m = v_t.unwrap().t().dot(&ss.dot(&u.unwrap().t()));
        println!("in {}s", now.elapsed().as_secs());

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
        m1.reset();
        m2.reset();
        t_xyz.y = 1e-6;
        m1.update(t_xyz, r_xyz, 1);
        r_xyz.x = 1e-6;
        m2.update(t_xyz, r_xyz, 2);
        gs.wavefront.reset();
        gs.reset_rays();
        m2.blocking(&mut gs.rays);
        m1.trace(&mut gs.rays);
        // gs.rays.gmt_truss_onaxis();
        m2.trace(&mut gs.rays);
        gs.rays.to_sphere1(-5.830, 2.197173);
        gs.opd2phase();
        wfs.propagate(&mut gs);
        wfs.process();
        let mut c: Vec<f32> = vec![0.0; n_c];
        dev2host(c.as_mut_ptr(), wfs.data_proc.d__c, n_c as i32);
        wfs.reset();

        let _s = Array::from_shape_vec((n_c, 1), c).unwrap();
        let _c_e = __m.dot(&_s);
        let _c_e2 = _c_e.into_shape((14, 6)).unwrap();
        println!("{:.2}", _c_e2 * 1e6);

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

            m1.reset();
            m2.reset();
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

            //println!("{:.2}", com[(0,0)] * 1e6);
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
            let mut id = 0;

            // GMT initial conditions
            let mut com = Array2::<f32>::zeros((14, 6));
            com[(0, 0)] = 1e-5;
            com[(8, 4)] = 1e-6;
            com[(2, 3)] = -1e-6;
            com[(11, 1)] = -1e-5;

            for k_step in 0..30 {
                println!("Step #{}", k_step);

                for sid in 1..8 {
                    //print!("{}", sid);
                    id = sid - 1;
                    t.x = com[[id, 0]] as f64;
                    t.y = com[[id, 1]] as f64;
                    t.z = com[[id, 2]] as f64;
                    r.x = com[[id, 3]] as f64;
                    r.y = com[[id, 4]] as f64;
                    r.z = com[[id, 5]] as f64;
                    m1.update(t, r, sid as i32);
                    id += 7;
                    t.x = com[[id, 0]] as f64;
                    t.y = com[[id, 1]] as f64;
                    t.z = com[[id, 2]] as f64;
                    r.x = com[[id, 3]] as f64;
                    r.y = com[[id, 4]] as f64;
                    r.z = com[[id, 5]] as f64;
                    m2.update(t, r, sid as i32)
                }
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
                println!(
                    "WFE RMS: {:.3} nm; PSSn: {}",
                    src_wfe_rms * 1e9,
                    src_pssn.eval()
                );

                gs.reset_rays();
                m2.blocking(&mut gs.rays);
                m1.trace(&mut gs.rays);
                // gs.rays.gmt_truss_onaxis();
                m2.trace(&mut gs.rays);
                gs.rays.to_sphere1(-5.830, 2.197173);
                gs.opd2phase();
                wfs.propagate(&mut gs);
                wfs.process();
                let mut c: Vec<f32> = vec![0.0; n_c];
                dev2host(c.as_mut_ptr(), wfs.data_proc.d__c, n_c as i32);
                wfs.reset();
                let _s = Array::from_shape_vec((n_c, 1), c).unwrap();
                let _c_e = __m.dot(&_s);
                let _c_e2 = _c_e.into_shape((14,6)).unwrap();

                //println!("{:.2}", c_e2 * 1e6);
                let q = gain * _c_e2;
//                let q_ = ArrayView2::from(&q);
//                com += &q.to_owned().view() ;
                com -= &q.view();
            }
        }

        m1_modes.cleanup();
        m1.cleanup();
        m2_modes.cleanup();
        m2.cleanup();
        gs.cleanup();
        wfs.cleanup();
    }
}
