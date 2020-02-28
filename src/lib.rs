#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use std::f32;
use std::time::{Instant};
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array2};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};

pub mod bindings;
pub mod ceo;

pub struct OpticalPathToSH48 {
    pub gmt: ceo::Gmt,
    pub gs: ceo::Source,
    pub sensor: ceo::GeometricShackHartmann,
}
impl OpticalPathToSH48 {
    pub fn new() -> OpticalPathToSH48 {
        let n_side_lenslet = 48;
        let n_px_lenslet = 16;
        //let n_px = n_side_lenslet*16 + 1;
        let pupil_size = 25.5;
        let n_gs = 3;
        let m1_n_mode = 27;
        let d = pupil_size / n_side_lenslet as f64;
        OpticalPathToSH48 {
            gmt: ceo::Gmt::new(m1_n_mode, None),
            gs: ceo::Source::empty(),
            sensor: ceo::Geometric_ShackHartmann::new(n_gs, n_side_lenslet, n_px_lenslet, d),
        }
    }
    pub fn build(&mut self) -> &mut Self {
        let z = 6. * f32::consts::PI / 180. / 60.;
        let a = 2. * f32::consts::PI / 3.;
        self.gmt.build();
        self.sensor.build();
        self.gs = self.sensor.new_guide_stars();
        self.gs.build(
            "V",
            vec![z, z, z],
            vec![0.0 * a, a, 2.0 * a],
            vec![0.0, 0.0, 0.0],
        );
        self.gs.through(&mut self.gmt);
        self.sensor.calibrate(&mut self.gs, 0.0).unwrap();
        self
    }
    pub fn propagate_src(&mut self) {
        self.gs.through(&mut self.gmt).through(&mut self.sensor);
    }
    pub fn calibrate(&mut self) -> Array2<f32> {
        let now = Instant::now();
        let n_rbm: usize = 84;
        let mut c_c_p: Vec<f32> = Vec::new();
        let mut c_c_m: Vec<f32> = Vec::new();
        let n_c: usize = (self.sensor.n_side_lenslet.pow(2) as usize) * 2 * self.sensor.n_sensor as usize;
        self.gmt.reset();
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
                            self.gmt.reset();
                            let mut t_xyz = vec![0.0; 3];
                            let mut r_xyz = vec![0.0; 3];
                            if tr == 1 {
                                t_xyz[a] = 1e-6;
                            }
                            if tr == 2 {
                                r_xyz[a] = 1e-6;
                            }
                            if mid == 1 {
                                self.gmt.set_m1_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            if mid == 2 {
                                self.gmt.set_m2_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            self.propagate_src();
                            self.sensor.process();
                            c_c_p.append(&mut self.sensor.centroids);
                        }
                    }
                }
            }
        }

        let pb_bm = ProgressBar::new(self.gmt.m1_n_mode as u64);
        pb_bm.set_style(ProgressStyle::default_bar().template("{bar:40.cyan/blue} {msg}"));
        if self.gmt.m1_n_mode > 0 {
            //println!("Bending modes ...");
            let mut a: Vec<f64> = vec![0.0; 7 * self.gmt.m1_n_mode as usize];
            self.gmt.reset();
            for sid in 0..7 {
                let s_ = String::from("S") + &(sid + 1).to_string() + "BM";
                pb_bm.set_message(&s_);
                pb_bm.reset();
                for k_a in 0..self.gmt.m1_n_mode {
                    pb_bm.inc(1);
                    let k = k_a + sid * self.gmt.m1_n_mode;
                    a[k as usize] = 1e-6;
                    self.gmt.set_m1_modes(&mut a);
                    self.propagate_src();
                    self.sensor.process();
                    c_c_p.append(&mut self.sensor.centroids);
                    a[k as usize] = 0.0;
                    self.gmt.set_m1_modes(&mut a);
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
                            self.gmt.reset();
                            let mut t_xyz = vec![0.0; 3];
                            let mut r_xyz = vec![0.0; 3];
                            if tr == 1 {
                                t_xyz[a] = -1e-6;
                            }
                            if tr == 2 {
                                r_xyz[a] = -1e-6;
                            }
                            if mid == 1 {
                                self.gmt.set_m1_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            if mid == 2 {
                                self.gmt.set_m2_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            self.propagate_src();
                            self.sensor.process();
                            //let sum2: f32 = c.iter().sum();
                            //println!("Centroids sum: {:.16}", sum2);
                            //println!("Diff. centroids sum: {:.14}", sum2 - sum1);
                            c_c_m.append(&mut self.sensor.centroids);
                        }
                    }
                }
            }
        }

        if self.gmt.m1_n_mode > 0 {
            //println!("Bending modes ...");
            let mut a: Vec<f64> = vec![0.0; 7 * self.gmt.m1_n_mode as usize];
            self.gmt.reset();
            for sid in 0..7 {
                let s_ = String::from("S") + &(sid + 1).to_string() + "BM";
                pb_bm.set_message(&s_);
                pb_bm.reset();
                for k_a in 0..self.gmt.m1_n_mode {
                    pb_bm.inc(1);
                    let k = k_a + sid * self.gmt.m1_n_mode;
                    a[k as usize] = -1e-6;
                    self.gmt.set_m1_modes(&mut a);
                    self.propagate_src();
                    self.sensor.process();
                    c_c_m.append(&mut self.sensor.centroids);
                    a[k as usize] = 0.0;
                    self.gmt.set_m1_modes(&mut a);
                }
            }
        }
        pb.finish();
        pb_bm.finish();

        //println!("WFS centroids #     : {}", self.sensor.n_centroids);
        //println!("WFS centroids length: {}", self.sensor.centroids.len());
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
        let n_sv = n_rbm + self.gmt.m1_n_mode as usize * 7;
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
        print!("Computing the pseudo-inverse");
        let now = Instant::now();
        let __m: Array2<f32> = _vt.t().dot(&l_sv.dot(&_u.t()));
        println!(" in {}ms", now.elapsed().as_millis());
        __m
    }
}
impl Drop for OpticalPathToSH48 {
    fn drop(&mut self) {
        drop(&mut self.gmt);
        drop(&mut self.gs);
        drop(&mut self.sensor);
    }
}
