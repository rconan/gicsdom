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
    pub fn calibrate_karhunen_loeve(
        &mut self,
        n_kl: usize,
        first_kl: Option<usize>,
        stroke: Option<f64>,
    ) -> Vec<Vec<f64>> {
        let mut gmt = ceo!(element::GMT, set_m2_n_mode = [n_kl]);
        let mut kl: Vec<Vec<f32>> = vec![];
        let first_kl = first_kl.unwrap_or(0);
        let stroke = stroke.unwrap_or(1e-6);
        for s in 0..7 {
            for k in first_kl..n_kl {
                gmt.set_m2_modes_ij(s, k, stroke);
                self.mmse_star.through(&mut gmt).xpupil();
                let b_push: Vec<f32> = self.mmse_star.phase_as_ptr().into();
                gmt.set_m2_modes_ij(s, k, -stroke);
                self.mmse_star.through(&mut gmt).xpupil();
                let b_pull: Vec<f32> = self.mmse_star.phase_as_ptr().into();
                kl.push(
                    b_push
                        .iter()
                        .zip(b_pull.iter())
                        .map(|x| 0.5 * (x.0 - x.1) / stroke as f32)
                        .collect::<Vec<f32>>(),
                );
            }
        }
        gmt.reset();
        let kl_norm = kl
            .iter()
            .map(|x| x.iter().map(|y| (y * y) as f64).sum::<f64>())
            .collect::<Vec<f64>>();
        kl.iter()
            .zip(kl_norm.into_iter())
            .map(|x| x.0.iter().map(|&y| y as f64 / x.1).collect::<Vec<f64>>())
            .collect::<_>()
    }
    pub fn get_karhunen_loeve_coefficients(
        &mut self,
        kln: &Vec<Vec<f64>>,
        stroke: Option<f64>,
    ) -> Vec<f64> {
        let stroke = stroke.unwrap_or(1f64);
        kln.iter()
            .map(|x| {
                x.iter()
                    .zip(Vec::<f32>::from(self.phase_as_ptr()).into_iter())
                    .map(|y| stroke * y.0 * y.1 as f64)
                    .sum::<f64>()
            })
            .collect::<_>()
    }
}
impl Drop for LinearMinimumMeanSquareError {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
