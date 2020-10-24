use super::ceo_bindings::LMMSE;
use super::{ceo, element, Atmosphere, Cu, GeometricShackHartmann as WFS, Mask, Source, CEO};
use std::ffi::CString;

pub struct LinearMinimumMeanSquareError {
    _c_: LMMSE,
    atm: Atmosphere,
    guide_star: Source,
    mmse_star: Source,
    fov_diameter: Option<f64>,
    pupil_mask: Mask,
}
impl CEO<element::LMMSE> {
    pub fn new() -> CEO<element::LMMSE> {
        CEO {
            args: element::LMMSE::default(),
        }
    }
    pub fn set_atmosphere(mut self, atm: &CEO<element::ATMOSPHERE>) -> Self {
        self.args.atm = atm.clone();
        self
    }
    pub fn set_guide_star(mut self, src: &CEO<element::SOURCE>) -> Self {
        self.args.guide_star = src.clone();
        self
    }
    pub fn set_mmse_star(mut self, src: &CEO<element::SOURCE>) -> Self {
        self.args.mmse_star = src.clone();
        self.args.fov_diameter = None;
        self
    }
    pub fn set_fov_diameter(mut self, fov_diameter: f64) -> Self {
        self.args.fov_diameter = Some(fov_diameter);
        self
    }
    pub fn set_n_side_lenslet(mut self, n_side_lenslet: usize) -> Self {
        self.args.n_side_lenslet = n_side_lenslet;
        self
    }
    pub fn build(self) -> LinearMinimumMeanSquareError {
        let mut gmt = ceo!(element::GMT);
        let mut mmse_star = self.args.mmse_star.build();
        mmse_star.through(&mut gmt).xpupil();
        let mut pupil_mask = Mask::new();
        let n_actuator = self.args.n_side_lenslet + 1;
        pupil_mask
            .build(n_actuator * n_actuator)
            .filter(&mut mmse_star.amplitude().into());
        let mut lmmse = LinearMinimumMeanSquareError {
            _c_: unsafe { std::mem::zeroed() },
            atm: self.args.atm.clone().build(),
            guide_star: self.args.guide_star.build(),
            mmse_star,
            fov_diameter: self.args.fov_diameter,
            pupil_mask,
        };
        let solver_id = CString::new(self.args.solver_id.clone().into_bytes()).unwrap();
        let d = self.args.guide_star.args.pupil_size / self.args.n_side_lenslet as f64;
        match lmmse.fov_diameter {
            Some(fov) => unsafe {
                lmmse._c_.setup3(
                    lmmse.atm.as_raw_mut_ptr(),
                    lmmse.guide_star.as_raw_mut_ptr(),
                    d as f32,
                    self.args.n_side_lenslet as i32,
                    lmmse.pupil_mask.as_raw_mut_ptr(),
                    solver_id.into_raw(),
                    self.args.wavefront_osf as i32,
                    0.5 * fov as f32,
                )
            },
            None => unsafe {
                lmmse._c_.setup2(
                    lmmse.atm.as_raw_mut_ptr(),
                    lmmse.guide_star.as_raw_mut_ptr(),
                    lmmse.mmse_star.as_raw_mut_ptr(),
                    d as f32,
                    self.args.n_side_lenslet as i32,
                    lmmse.pupil_mask.as_raw_mut_ptr(),
                    solver_id.into_raw(),
                    self.args.wavefront_osf as i32,
                )
            },
        }
        lmmse
    }
}
impl LinearMinimumMeanSquareError {
    pub fn get_wavefront_estimate(&mut self, wfs: &mut WFS) -> &mut Self {
        unsafe {
            self._c_.estimation(&wfs.as_raw_mut_ptr().data_proc);
        }
        self
    }
    pub fn phase_as_ptr(&mut self) -> Cu<f32> {
        //println!("PS_E_N_PX: {}",self._c_.PS_E_N_PX);
        let mut phase: Cu<f32> = Cu::vector(self._c_.PS_E_N_PX as usize);
        phase.from_ptr(self._c_.d__phase_est);
        phase
    }
    pub fn set_n_iteration(&mut self, n_iteration: usize) {
        self._c_.iSolve.N_ITERATION = n_iteration as i32;
    }
    pub fn get_n_iteration(&mut self) -> usize {
        self._c_.iSolve.N_ITERATION as usize
    }
}
