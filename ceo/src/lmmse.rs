use super::ceo_bindings::LMMSE;
use super::{element, Atmosphere, Mask, Source, CEO};
use std::ffi::CString;

pub struct LinearMinimumMeanSquareError {
    _c_: LMMSE,
    atm: Atmosphere,
    guide_star: Source,
    mmse_star: Option<Source>,
    fov_diameter: Option<f64>,
    pupil_mask: Mask,
}
impl CEO<element::LMMSE> {
    pub fn new() -> CEO<element::LMMSE> {
        CEO {
            args: element::LMMSE::default(),
        }
    }
    pub fn set_atm(mut self, atm: CEO<element::ATMOSPHERE>) -> Self {
        self.args.atm = atm;
        self
    }
    pub fn set_guide_star(mut self, src: CEO<element::SOURCE>) -> Self {
        self.args.guide_star = src;
        self
    }
    pub fn set_mmse_star(mut self, src: CEO<element::SOURCE>) -> Self {
        self.args.mmse_star = Some(src);
        self.args.fov_diameter = None;
        self
    }
    pub fn set_fov_diameter(mut self, fov_diameter: f64) -> Self {
        self.args.fov_diameter = Some(fov_diameter);
        self.args.mmse_star = None;
        self
    }
    pub fn set_n_side_lenslet(mut self, n_side_lenslet: usize) -> Self {
        self.args.n_side_lenslet = n_side_lenslet;
        self
    }
    pub fn set_pupil_mask(mut self, pupil_mask: Mask) -> Self {
        self.args.pupil_mask = pupil_mask;
        self
    }
    pub fn build(self) -> LinearMinimumMeanSquareError {
        let mut lmmse = LinearMinimumMeanSquareError {
            _c_: unsafe { std::mem::zeroed() },
            atm: self.args.atm.clone().build(),
            guide_star: self.args.guide_star.build(),
            mmse_star: self.args.mmse_star.as_ref().map(|x| x.build()),
            fov_diameter: self.args.fov_diameter,
            pupil_mask: self.args.pupil_mask,
        };
        let solver_id = CString::new(self.args.solver_id.clone().into_bytes()).unwrap();
        match lmmse.fov_diameter {
            Some(fov) => unsafe {
                lmmse._c_.setup3(
                    lmmse.atm.as_raw_mut_ptr(),
                    lmmse.guide_star.as_raw_mut_ptr(),
                    self.args.guide_star.args.pupil_sampling as f32,
                    self.args.n_side_lenslet as i32,
                    lmmse.pupil_mask.as_raw_mut_ptr(),
                    solver_id.into_raw(),
                    self.args.wavefront_osf as i32,
                    0.5 * fov as f32,
                )
            },
            None => {
                if let Some(ref mut src) = lmmse.mmse_star {
                    unsafe {
                        lmmse._c_.setup2(
                            lmmse.atm.as_raw_mut_ptr(),
                            lmmse.guide_star.as_raw_mut_ptr(),
                            src.as_raw_mut_ptr(),
                            self.args.guide_star.args.pupil_sampling as f32,
                            self.args.n_side_lenslet as i32,
                            lmmse.pupil_mask.as_raw_mut_ptr(),
                            solver_id.into_raw(),
                            self.args.wavefront_osf as i32,
                        )
                    }
                }
            }
        }
        lmmse
    }
}
