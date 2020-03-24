use std::ffi::CString;
use std::{f32, mem};

use super::ceo_bindings::{source, vector};

pub struct Source {
    pub _c_: source,
    size: i32,
    pub pupil_size: f64,
    pub pupil_sampling: i32,
    pub _wfe_rms: Vec<f32>,
}
impl Source {
    pub fn empty() -> Source {
        Source {
            _c_: unsafe { mem::zeroed() },
            size: 0,
            pupil_size: 0.0,
            pupil_sampling: 0,
            _wfe_rms: vec![],
        }
    }
    pub fn new(size: i32, pupil_size: f64, pupil_sampling: i32) -> Source {
        Source {
            _c_: unsafe { mem::zeroed() },
            size,
            pupil_size,
            pupil_sampling,
            _wfe_rms: vec![0.0; size as usize],
        }
    }
    pub fn from(args: (i32, f64, i32)) -> Source {
        Source::new(args.0, args.1, args.2)
    }
    pub fn build(
        &mut self,
        band: &str,
        mut zenith: Vec<f32>,
        mut azimuth: Vec<f32>,
        mut magnitude: Vec<f32>,
    ) -> &mut Self {
        assert_eq!(zenith.len(), azimuth.len());
        assert_eq!(zenith.len(), magnitude.len());
        let band = CString::new(band).unwrap();
        unsafe {
            let origin = vector {
                x: 0.0,
                y: 0.0,
                z: 25.0,
            };
            self._c_.setup7(
                band.into_raw(),
                magnitude.as_mut_ptr(),
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                f32::INFINITY,
                self.size,
                self.pupil_size,
                self.pupil_sampling,
                origin,
            );
        }
        self
    }
    pub fn set_fwhm(&mut self, value: f64) {
        self._c_.fwhm = value as f32;
    }
    pub fn rotate_rays(&mut self, angle: f64) {
        self._c_.rays.rot_angle = angle;
    }
    pub fn xpupil(&mut self) -> &mut Self {
        unsafe {
            self._c_.wavefront.reset();
            self._c_.opd2phase();
        }
        self
    }
    pub fn wfe_rms(&mut self) -> Vec<f32> {
        unsafe {
            self._c_.wavefront.rms(self._wfe_rms.as_mut_ptr());
        }
        self._wfe_rms.clone()
    }
    pub fn wfe_rms_10e(&mut self, exp: i32) -> Vec<f32> {
        unsafe {
            self._c_.wavefront.rms(self._wfe_rms.as_mut_ptr());
        }
        self._wfe_rms
            .iter()
            .map(|x| x * 10_f32.powi(-exp))
            .collect()
    }
    pub fn segments_gradients(&mut self) -> Vec<Vec<f32>> {
        let mut sxy: Vec<Vec<f32>> = vec![vec![0.;7 * self.size as usize]; 2];
        unsafe {
            self._c_.wavefront.segments_gradient_averageFast(
                sxy[0].as_mut_ptr(),
                sxy[1].as_mut_ptr(),
                self._c_.rays.L as f32,
                self._c_.rays.d__piston_mask,
            );
        }
        sxy
    }
    pub fn reset(&mut self) {
        unsafe {
            self._c_.wavefront.reset();
            self._c_.reset_rays();
        }
    }
    pub fn update(&mut self, mut zenith: Vec<f64>, mut azimuth: Vec<f64>) {
        unsafe {
            self._c_.update_directions(
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                zenith.len() as i32,
            );
        }
    }
}
impl Drop for Source {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
pub trait Propagation {
    fn propagate(&mut self, src: &mut Source) -> &mut Self;
    fn time_propagate(&mut self, secs: f64, src: &mut Source) -> &mut Self;
}
impl Source {
    pub fn through<T: Propagation>(&mut self, system: &mut T) -> &mut Self {
        system.propagate(self);
        self
    }
}
