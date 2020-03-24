use ndarray::Array2;
use std::ffi::CString;
use std::{f32, mem};

use super::ceo_bindings::{bundle, gmt_m1, gmt_m2, modes, vector, zernikeS};
use super::Propagation;
use super::Source;

#[derive(Clone, Debug)]
pub struct GmtState {
    pub rbm: Array2<f32>,
    pub bm: Array2<f32>,

}
pub struct Gmt {
    _c_m1_modes: modes,
    _c_m2_modes: zernikeS,
    _c_m1: gmt_m1,
    _c_m2: gmt_m2,
    pub m1_n_mode: u64,
    pub m2_n_mode: u64,
    a: Vec<f64>,
}
impl Gmt {
    pub fn test(m1_n_mode: u64, m2_n_mode: Option<u64>) -> Gmt {
        let mode_type = CString::new("bending modes").unwrap();
        let mut this: Gmt = Gmt {
            _c_m1_modes: unsafe { mem::zeroed() },
            _c_m2_modes: unsafe { mem::zeroed() },
            _c_m1: unsafe { mem::zeroed() },
            _c_m2: unsafe { mem::zeroed() },
            m1_n_mode,
            m2_n_mode: match m2_n_mode {
                Some(m2_n_mode) => m2_n_mode,
                None => 0,
            },
            a: vec![0.0],
        };
        if this.m2_n_mode > 0 {
            this.a = vec![0.0; this.m2_n_mode as usize];
        }
        unsafe {
            this._c_m1_modes
                .setup(mode_type.into_raw(), 7, this.m1_n_mode as i32);
            this._c_m1.setup1(&mut this._c_m1_modes);
            this._c_m2_modes
                .setup1(this.m2_n_mode as i32, this.a.as_mut_ptr(), 7);
            this._c_m2.setup2(&mut this._c_m2_modes);
        }
        this
    }
    pub fn new() -> Gmt {
        Gmt {
            _c_m1_modes: unsafe { mem::zeroed() },
            _c_m2_modes: unsafe { mem::zeroed() },
            _c_m1: unsafe { mem::zeroed() },
            _c_m2: unsafe { mem::zeroed() },
            m1_n_mode: 0,
            m2_n_mode: 0,
            a: vec![0.],
        }
    }
    pub fn build(&mut self,m1_n_mode: u64, m2_n_mode: Option<u64>) -> &mut Gmt {
        let mode_type = CString::new("bending modes").unwrap();
        self.m1_n_mode = m1_n_mode;
        self.m2_n_mode = match m2_n_mode {
            Some(m2_n_mode) => m2_n_mode,
            None => 0,
        };
        if self.m2_n_mode > 0 {
            self.a = vec![0.0; self.m2_n_mode as usize];
        }
        unsafe {
            self._c_m1_modes
                .setup(mode_type.into_raw(), 7, self.m1_n_mode as i32);
            self._c_m1.setup1(&mut self._c_m1_modes);
            self._c_m2_modes
                .setup1(self.m2_n_mode as i32, self.a.as_mut_ptr(), 7);
            self._c_m2.setup2(&mut self._c_m2_modes);
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        let mut a: Vec<f64> = vec![0.0; 7 * self.m1_n_mode as usize];
        unsafe {
            self._c_m1.reset();
            self._c_m2.reset();
            self._c_m1_modes.update(a.as_mut_ptr());
        }
        self
    }
    pub fn set_m1_segment_state(&mut self, sid: i32, t_xyz: &Vec<f64>, r_xyz: &Vec<f64>) {
        let t_xyz = vector {
            x: t_xyz[0],
            y: t_xyz[1],
            z: t_xyz[2],
        };
        let r_xyz = vector {
            x: r_xyz[0],
            y: r_xyz[1],
            z: r_xyz[2],
        };
        unsafe {
            self._c_m1.update(t_xyz, r_xyz, sid);
        }
    }
    pub fn set_m2_segment_state(&mut self, sid: i32, t_xyz: &Vec<f64>, r_xyz: &Vec<f64>) {
        let t_xyz = vector {
            x: t_xyz[0],
            y: t_xyz[1],
            z: t_xyz[2],
        };
        let r_xyz = vector {
            x: r_xyz[0],
            y: r_xyz[1],
            z: r_xyz[2],
        };
        unsafe {
            self._c_m2.update(t_xyz, r_xyz, sid);
        }
    }
    pub fn set_m1_modes(&mut self, a: &mut Vec<f64>) {
        unsafe {
            self._c_m1_modes.update(a.as_mut_ptr());
        }
    }
    pub fn update(&mut self, gstate: &GmtState) {
        let mut t_xyz = vec![0.0; 3];
        let mut r_xyz = vec![0.0; 3];
        let mut a: Vec<f64> = vec![0.0; 7 * self.m1_n_mode as usize];
        let mut id = 0;

        for sid in 1..8 {
            //print!("{}", sid);••••••••••••
            id = sid - 1;
            t_xyz[0] = gstate.rbm[[id, 0]] as f64;
            t_xyz[1] = gstate.rbm[[id, 1]] as f64;
            t_xyz[2] = gstate.rbm[[id, 2]] as f64;
            r_xyz[0] = gstate.rbm[[id, 3]] as f64;
            r_xyz[1] = gstate.rbm[[id, 4]] as f64;
            r_xyz[2] = gstate.rbm[[id, 5]] as f64;
            self.set_m1_segment_state(sid as i32, &t_xyz, &r_xyz);
            if self.m1_n_mode > 0 {
                for k_bm in 0..self.m1_n_mode {
                    let idx = id * self.m1_n_mode as usize + k_bm as usize;
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
            self.set_m2_segment_state(sid as i32, &t_xyz, &r_xyz);
        }
        self.set_m1_modes(&mut a);
    }
}
impl Drop for Gmt {
    fn drop(&mut self) {
        unsafe {
            self._c_m1_modes.cleanup();
            self._c_m1.cleanup();
            self._c_m2_modes.cleanup();
            self._c_m2.cleanup();
        }
    }
}
impl Propagation for Gmt {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            src._c_.reset_rays();
            let rays: &mut bundle = &mut src._c_.rays;
            self._c_m2.blocking(rays);
            self._c_m1.trace(rays);
            //gs.rays.gmt_truss_onaxis();
            self._c_m2.trace(rays);
            rays.to_sphere1(-5.830, 2.197173);
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}
