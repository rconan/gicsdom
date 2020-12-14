//!
//! # CEO gmt wrapper
//!
//! Provides a structure `Gmt` that is a wrapper for [CEO](https://github.com/rconan/CEO) gmt C++ structure.
//! `Gmt` is instantiated and initialized with the `GMT` builder
//!
//! # Examples
//!
//! ```
//! use ceo::ceo;
//! // Creates a gmt instance with default parameters
//! let mut src = ceo!(GMT);
//! ```
//!
//! ```
//! use ceo::ceo;
//! // Creates a gmt instance with 27 M1 bending modes
//! let mut gmt = ceo!(GMT, set_m1_n_mode = [27]);
//! ```

use super::ceo_bindings::{bundle, gmt_m1, gmt_m2, modes, vector};
use super::{Builder, Propagation, Source};
use std::{ffi::{CString,CStr}, mem};

#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct Mirror {
    pub mode_type: String,
    pub n_mode: usize,
}
/// `Gmt` builder
///
/// Default properties:
///  - M1:
///    - mode type : "bending modes"
///    - \# mode    : 0
///  - M2:
///    - mode type : "Karhunen-Loeve"
///    - \# mode    : 0
///
/// # Examples
///
/// ```
/// use ceo::{Builder, GMT};
/// let mut src = GMT::new().build();
/// ```
///
/// ```
/// use ceo::{Builder, GMT};
/// let mut gmt = GMT::new().set_m1_n_mode(27).build();
/// ```
#[derive(Debug, Clone)]
pub struct GMT {
    pub m1: Mirror,
    pub m2: Mirror,
}
impl Default for GMT {
    fn default() -> Self {
        GMT {
            m1: Mirror {
                mode_type: "bending modes".into(),
                n_mode: 0,
            },
            m2: Mirror {
                mode_type: "Karhunen-Loeve".into(),
                n_mode: 0,
            },
        }
    }
}
impl GMT {
    pub fn set_m1(self, mode_type: &str, n_mode: usize) -> Self {
        Self {
            m1: Mirror {
                mode_type: mode_type.into(),
                n_mode,
            },
            ..self
        }
    }
    pub fn set_m1_n_mode(self, n_mode: usize) -> Self {
        Self {
            m1: Mirror { n_mode, ..self.m1 },
            ..self
        }
    }
    pub fn set_m2(self, mode_type: &str, n_mode: usize) -> Self {
        Self {
            m2: Mirror {
                mode_type: mode_type.into(),
                n_mode,
            },
            ..self
        }
    }
    pub fn set_m2_n_mode(self, n_mode: usize) -> Self {
        Self {
            m2: Mirror { n_mode, ..self.m2 },
            ..self
        }
    }
}
impl Builder for GMT {
    type Component = Gmt;
    fn build(self) -> Gmt {
        let mut gmt = Gmt {
            _c_m1_modes: unsafe { mem::zeroed() },
            _c_m2_modes: unsafe { mem::zeroed() },
            _c_m1: unsafe { mem::zeroed() },
            _c_m2: unsafe { mem::zeroed() },
            m1_n_mode: 0,
            m2_n_mode: 0,
            m2_max_n: 0,
            a1: vec![0.],
            a2: vec![0.],
        };
        let m1_mode_type = CString::new(self.m1.mode_type.into_bytes()).unwrap();
        gmt.m1_n_mode = self.m1.n_mode;
        gmt.a1 = vec![0.0; 7 * gmt.m1_n_mode as usize];
        unsafe {
            gmt._c_m1_modes
                .setup(m1_mode_type.into_raw(), 7, gmt.m1_n_mode as i32);
            gmt._c_m1.setup1(&mut gmt._c_m1_modes);
        }
        let m2_mode_type = CString::new(self.m2.mode_type.into_bytes()).unwrap();
        gmt.m2_n_mode = self.m2.n_mode;
        gmt.a2 = vec![0.0; 7 * gmt.m2_n_mode as usize];
        unsafe {
            gmt._c_m2_modes
                .setup(m2_mode_type.into_raw(), 7, gmt.m2_n_mode as i32);
            gmt._c_m2.setup1(&mut gmt._c_m2_modes);
        }
        gmt
    }
}
impl From<&Gmt> for GMT {
    fn from(gmt: &Gmt) -> Self {
        Self { m1: gmt.get_m1(), m2: gmt.get_m2()}
    }
}
/// gmt wrapper
pub struct Gmt {
    pub _c_m1_modes: modes,
    pub _c_m2_modes: modes,
    pub _c_m1: gmt_m1,
    pub _c_m2: gmt_m2,
    /// M1 number of bending modes per segment
    pub m1_n_mode: usize,
    /// M2 number of bending modes per segment
    pub m2_n_mode: usize,
    /// M2 largest Zernike radial order per segment
    pub m2_max_n: usize,
    pub a1: Vec<f64>,
    pub a2: Vec<f64>,
}
impl Gmt {
    /// Creates a new `Gmt`
    pub fn new() -> Gmt {
        Gmt {
            _c_m1_modes: unsafe { mem::zeroed() },
            _c_m2_modes: unsafe { mem::zeroed() },
            _c_m1: unsafe { mem::zeroed() },
            _c_m2: unsafe { mem::zeroed() },
            m1_n_mode: 0,
            m2_n_mode: 1,
            m2_max_n: 0,
            a1: vec![0.],
            a2: vec![0.],
        }
    }
    /// Sets the `Gmt` parameters:
    ///
    /// * `m1_n_mode` - the number of of modes on each M1 segment
    /// * `m2_max_n` - M2 largest Zernike radial order per segment
    pub fn build(&mut self, m1_n_mode: usize, m2_max_n: Option<usize>) -> &mut Gmt {
        let mode_type = CString::new("bending modes").unwrap();
        self.m1_n_mode = m1_n_mode;
        self.m2_max_n = match m2_max_n {
            Some(m2_max_n) => m2_max_n,
            None => 0,
        };
        self.m2_n_mode = (self.m2_max_n + 1) * (self.m2_max_n + 2) / 2;
        if self.m1_n_mode > 0 {
            self.a1 = vec![0.0; 7 * self.m1_n_mode as usize];
        }
        self.a2 = vec![0.0; 7 * self.m2_n_mode as usize];
        unsafe {
            self._c_m1_modes
                .setup(mode_type.into_raw(), 7, self.m1_n_mode as i32);
            self._c_m1.setup1(&mut self._c_m1_modes);
            self._c_m2_modes.setup(
                CString::new("Karhunen-Loeve").unwrap().into_raw(),
                7,
                self.m2_n_mode as i32,
            );
            self._c_m2.setup1(&mut self._c_m2_modes);
        }
        self
    }
    /// Sets the `Gmt` M1 parameters:
    ///
    /// * `mode_type` - the type of modes: "bending modes", "KarhunenLoeve", ...
    /// * `m1_n_mode` - the number of modes on each M1 segment
    pub fn build_m1(&mut self, mode_type: &str, m1_n_mode: usize) -> &mut Self {
        let m1_mode_type = CString::new(mode_type).unwrap();
        self.m1_n_mode = m1_n_mode;
        if self.m1_n_mode > 0 {
            self.a1 = vec![0.0; 7 * self.m1_n_mode as usize];
        }
        unsafe {
            self._c_m1_modes
                .setup(m1_mode_type.into_raw(), 7, self.m1_n_mode as i32);
            self._c_m1.setup1(&mut self._c_m1_modes);
        }
        self
    }
    pub fn from_m2_modes(&mut self, mode_type: &str, m2_n_mode: usize) -> &mut Self {
        let m2_mode_type = CString::new(mode_type).unwrap();
        self.m2_n_mode = m2_n_mode;
        if self.m2_n_mode > 0 {
            self.a1 = vec![0.0; 7 * self.m2_n_mode as usize];
        }
        unsafe {
            self._c_m2_modes
                .setup(m2_mode_type.into_raw(), 7, self.m2_n_mode as i32);
            self._c_m2.setup1(&mut self._c_m2_modes);
        }
        self
    }
    /// Sets the `Gmt` M2 parameters:
    ///
    /// * `m2_n_mode` - the number of Karhunen-Loeve modes on each M2 segment
    pub fn build_m2(&mut self, m2_n_mode: Option<usize>) -> &mut Self {
        let m2_mode_type = CString::new("Karhunen-Loeve").unwrap();
        self.m2_n_mode = match m2_n_mode {
            Some(m2_n_mode) => m2_n_mode,
            None => 0,
        };
        self.a2 = vec![0.0; 7 * self.m2_n_mode as usize];
        unsafe {
            self._c_m2_modes
                .setup(m2_mode_type.into_raw(), 7, self.m2_n_mode as i32);
            self._c_m2.setup1(&mut self._c_m2_modes);
        }
        self
    }
    /// Returns `Gmt` M1 mode type
    pub fn get_m1_mode_type(&self) -> String {
        unsafe {
            String::from(
                CStr::from_ptr(self._c_m1_modes.filename.as_ptr())
                    .to_str()
                    .expect("CStr::to_str failed"),
            )
        }
    }
    /// Returns `Gmt` M1 properties
    pub fn get_m1(&self) -> Mirror {
        Mirror {
            mode_type: self.get_m1_mode_type(),
            n_mode: self.m1_n_mode,
        }
    }
    /// Returns `Gmt` M2 properties
    pub fn get_m2(&self) -> Mirror {
        Mirror {
            mode_type: self.get_m2_mode_type(),
            n_mode: self.m2_n_mode,
        }
    }
    /// Returns `Gmt` M2 mode type
    pub fn get_m2_mode_type(&self) -> String {
        unsafe {
            String::from(
                CStr::from_ptr(self._c_m2_modes.filename.as_ptr())
                    .to_str()
                    .expect("CStr::to_str failed"),
            )
        }
    }
    /// Resets M1 and M2 to their aligned states
    pub fn reset(&mut self) -> &mut Self {
        let mut a1: Vec<f64> = vec![0.0; 7 * self.m1_n_mode as usize];
        let mut a2: Vec<f64> = vec![0.0; 7 * self.m2_n_mode as usize];
        unsafe {
            self._c_m1.reset();
            self._c_m2.reset();
            self._c_m1_modes.update(a1.as_mut_ptr());
            self._c_m2_modes.update(a2.as_mut_ptr());
        }
        self
    }
    /// Keeps only the M1 segment specified in the vector `sid`
    ///
    /// * `sid` - vector of segment ID numbers in the range [1,7]
    pub fn keep(&mut self, sid: &mut Vec<i32>) -> &mut Self {
        unsafe {
            self._c_m1.keep(sid.as_mut_ptr(), sid.len() as i32);
            self._c_m2.keep(sid.as_mut_ptr(), sid.len() as i32);
        }
        self
    }
    /// Sets M1 segment rigid body motion with:
    ///
    /// * `sid` - the segment ID number in the range [1,7]
    /// * `t_xyz` - the 3 translations Tx, Ty and Tz
    /// * `r_xyz` - the 3 rotations Rx, Ry and Rz
    pub fn set_m1_segment_state(&mut self, sid: i32, t_xyz: &[f64], r_xyz: &[f64]) {
        assert!(sid > 0 && sid < 8, "Segment ID must be in the range [1,7]!");
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
    /// Sets M2 segment rigid body motion with:
    ///
    /// * `sid` - the segment ID number in the range [1,7]
    /// * `t_xyz` - the 3 translations Tx, Ty and Tz
    /// * `r_xyz` - the 3 rotations Rx, Ry and Rz
    pub fn set_m2_segment_state(&mut self, sid: i32, t_xyz: &[f64], r_xyz: &[f64]) {
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
    /// Sets M1 modal coefficients
    pub fn set_m1_modes(&mut self, a: &mut Vec<f64>) {
        unsafe {
            self._c_m1_modes.update(a.as_mut_ptr());
        }
    }
    pub fn set_m1_modes_ij(&mut self, i: usize, j: usize, value: f64) {
        let mut a = vec![0f64; 7 * self.m1_n_mode];
        a[i * self.m1_n_mode + j] = value;
        unsafe {
            self._c_m1_modes.update(a.as_mut_ptr());
        }
    }
    /// Sets M2 modal coefficients
    pub fn set_m2_modes(&mut self, a: &mut Vec<f64>) {
        unsafe {
            self._c_m2_modes.update(a.as_mut_ptr());
        }
    }
    pub fn set_m2_modes_ij(&mut self, i: usize, j: usize, value: f64) {
        let mut a = vec![0f64; 7 * self.m2_n_mode];
        a[i * self.m2_n_mode + j] = value;
        unsafe {
            self._c_m2_modes.update(a.as_mut_ptr());
        }
    }
    /// Updates M1 and M1 rigid body motion and M1 model coefficients
    pub fn update(
        &mut self,
        m1_rbm: Option<&Vec<Vec<f64>>>,
        m2_rbm: Option<&Vec<Vec<f64>>>,
        m1_mode: Option<&Vec<Vec<f64>>>,
        m2_mode: Option<&Vec<Vec<f64>>>,
    ) {
        if let Some(m1_rbm) = m1_rbm {
            for (k, rbm) in m1_rbm.iter().enumerate() {
                self.set_m1_segment_state((k + 1) as i32, &rbm[..3], &rbm[3..]);
            }
        }
        if let Some(m2_rbm) = m2_rbm {
            for (k, rbm) in m2_rbm.iter().enumerate() {
                self.set_m2_segment_state((k + 1) as i32, &rbm[..3], &rbm[3..]);
            }
        }
        if let Some(m1_mode) = m1_mode {
            let mut m = m1_mode.clone().into_iter().flatten().collect::<Vec<f64>>();
            self.set_m1_modes(&mut m);
        }
        if let Some(m2_mode) = m2_mode {
            let mut m = m2_mode.clone().into_iter().flatten().collect::<Vec<f64>>();
            self.set_m2_modes(&mut m);
        }
    }
    pub fn update42(
        &mut self,
        m1_rbm: Option<&Vec<f64>>,
        m2_rbm: Option<&Vec<f64>>,
        m1_mode: Option<&Vec<f64>>,
        m2_mode: Option<&Vec<f64>>,
    ) {
        if let Some(m1_rbm) = m1_rbm {
            for (k, rbm) in m1_rbm.chunks(6).enumerate() {
                self.set_m1_segment_state((k + 1) as i32, &rbm[..3], &rbm[3..]);
            }
        }
        if let Some(m2_rbm) = m2_rbm {
            for (k, rbm) in m2_rbm.chunks(6).enumerate() {
                self.set_m2_segment_state((k + 1) as i32, &rbm[..3], &rbm[3..]);
            }
        }
        if let Some(m1_mode) = m1_mode {
            let mut m = m1_mode.clone();
            self.set_m1_modes(&mut m);
        }
        if let Some(m2_mode) = m2_mode {
            let mut m = m2_mode.clone();
            self.set_m2_modes(&mut m);
        }
    }
    /*
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
    */
}
impl Drop for Gmt {
    /// Frees CEO memory before dropping `Gmt`
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
    /// Ray traces a `Source` through `Gmt`, ray tracing stops at the exit pupil
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            src.as_raw_mut_ptr().reset_rays();
            let rays: &mut bundle = &mut src.as_raw_mut_ptr().rays;
            self._c_m2.blocking(rays);
            self._c_m1.trace(rays);
            rays.gmt_truss_onaxis();
            rays.gmt_m2_baffle();
            self._c_m2.trace(rays);
            rays.to_sphere1(-5.830, 2.197173);
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Builder;

    #[test]
    fn gmt_new() {
        GMT::new().set_m1_n_mode(27).set_m2_n_mode(123).build();
    }

    #[test]
    fn gmt_new_with_macro() {
        crate::ceo!(GMT, set_m1_n_mode = [27], set_m2_n_mode = [123]);
    }

    #[test]
    fn gmt_optical_alignment() {
        use crate::SOURCE;
        let mut src = SOURCE::new().set_pupil_sampling(1001).build();
        let mut gmt = GMT::new().build();
        src.through(&mut gmt).xpupil();
        assert!(src.wfe_rms_10e(-9)[0] < 1.0);
    }

    #[test]
    fn gmt_m1_rx_optical_sensitity() {
        use crate::SOURCE;
        let mut src = SOURCE::new().set_pupil_sampling(1001).build();
        let mut gmt = GMT::new().build();
        let seg_tts0 = src.through(&mut gmt).xpupil().segments_gradients();
        let rt = vec![vec![0f64, 0f64, 0f64, 1e-6, 0f64, 0f64]; 7];
        gmt.update(Some(&rt), None, None, None);
        let seg_tts = src.through(&mut gmt).xpupil().segments_gradients();
        let mut delta: Vec<f32> = Vec::with_capacity(7);
        for k in 0..7 {
            delta
                .push(1e6 * (seg_tts[0][k] - seg_tts0[0][k]).hypot(seg_tts[1][k] - seg_tts0[1][k]));
        }
        assert!(delta.iter().all(|x| (x - 2.0).abs() < 1e-1));
    }

    #[test]
    fn gmt_m1_ry_optical_sensitity() {
        use crate::SOURCE;
        let mut src = SOURCE::new().set_pupil_sampling(1001).build();
        let mut gmt = GMT::new().build();
        let seg_tts0 = src.through(&mut gmt).xpupil().segments_gradients();
        let rt = vec![vec![0f64, 0f64, 0f64, 0f64, 1e-6, 0f64]; 7];
        gmt.update(Some(&rt), None, None, None);
        let seg_tts = src.through(&mut gmt).xpupil().segments_gradients();
        let mut delta: Vec<f32> = Vec::with_capacity(7);
        for k in 0..7 {
            delta
                .push(1e6 * (seg_tts[0][k] - seg_tts0[0][k]).hypot(seg_tts[1][k] - seg_tts0[1][k]));
        }
        assert!(delta.iter().all(|x| (x - 2.0).abs() < 1e-1));
    }

    #[test]
    fn gmt_m2_rx_optical_sensitity() {
        use crate::SOURCE;
        let mut src = SOURCE::new().set_pupil_sampling(1001).build();
        let mut gmt = GMT::new().build();
        let seg_tts0 = src.through(&mut gmt).xpupil().segments_gradients();
        let rt = vec![vec![0f64, 0f64, 0f64, 1e-6, 0f64, 0f64]; 7];
        gmt.update(None, Some(&rt), None, None);
        let seg_tts = src.through(&mut gmt).xpupil().segments_gradients();
        let mut delta: Vec<f32> = Vec::with_capacity(7);
        for k in 0..7 {
            delta
                .push(1e6 * (seg_tts[0][k] - seg_tts0[0][k]).hypot(seg_tts[1][k] - seg_tts0[1][k]));
        }
        assert!(delta.iter().all(|x| (x - 0.25).abs() < 1e-3));
    }

    #[test]
    fn gmt_m2_ry_optical_sensitity() {
        use crate::SOURCE;
        let mut src = SOURCE::new().set_pupil_sampling(1001).build();
        let mut gmt = GMT::new().build();
        let seg_tts0 = src.through(&mut gmt).xpupil().segments_gradients();
        let rt = vec![vec![0f64, 0f64, 0f64, 0f64, 1e-6, 0f64]; 7];
        gmt.update(None, Some(&rt), None, None);
        let seg_tts = src.through(&mut gmt).xpupil().segments_gradients();
        let mut delta: Vec<f32> = Vec::with_capacity(7);
        for k in 0..7 {
            delta
                .push(1e6 * (seg_tts[0][k] - seg_tts0[0][k]).hypot(seg_tts[1][k] - seg_tts0[1][k]));
        }
        assert!(delta.iter().all(|x| (x - 0.25).abs() < 1e-2));
    }

    /*
    #[test]
    fn gmt_lenslet_gradients() {
        let pupil_size = 25.5f64;
        let n_lenslet = 48i32;
        let lenslet_size = pupil_size / n_lenslet as f64;
        let mut gmt = Gmt::new();
        gmt.build(1, None);
        let mut src = Source::new(1, pupil_size, n_lenslet * 16 + 1);
        src.build("V", vec![0.0], vec![0.0], vec![0.0]);
        src.set_fwhm(4.0);
        let mut sensor = Imaging::new();
        sensor.build(1, n_lenslet as i32, 16, 2, 24, 3);
        let mut cog0 = Centroiding::new();
        cog0.build(n_lenslet as u32, None);

        src.through(&mut gmt).xpupil().through(&mut sensor);
        cog0.process(&sensor, None)
            .set_valid_lenslets(Some(0.9), None);
        src.lenslet_gradients(n_lenslet, lenslet_size, &mut cog0);
        let s0 = cog0.grab().valids(None);
        if s0.iter().any(|x| x.is_nan()) {
            let n = (n_lenslet * n_lenslet) as usize;
            for k in 0..n {
                if k % n_lenslet as usize == 0 {
                    println!("");
                }
                if cog0.centroids[k].is_nan() || cog0.centroids[k + n].is_nan() {
                    print!("X");
                } else {
                    print!("o");
                }
            }
        }
        let s0_any_nan = s0.iter().any(|x| x.is_nan());
        assert!(!s0_any_nan);

        let rt = vec![vec![0f64, 0f64, 0f64, 1e-6, 1e-6, 0f64]; 7];
        gmt.update(None, Some(&rt), None);

        let mut cog = Centroiding::new();
        cog.build(n_lenslet as u32, None);
        cog.set_valid_lenslets(None, Some(cog0.valid_lenslets.clone()));
        src.through(&mut gmt)
            .xpupil()
            .lenslet_gradients(n_lenslet, lenslet_size, &mut cog);

        let s = cog.grab().valids(None);
        let s_any_nan = s.iter().any(|x| x.is_nan());
        assert!(!s_any_nan);
    }
    */
}
