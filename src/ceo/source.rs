use std::ffi::CString;
use std::{f32, mem};

use super::ceo_bindings::{source, vector};
use super::{Centroiding, CuFloat};

/// A system that mutates `Source` arguments should implement the `Propagation` trait
pub trait Propagation {
    fn propagate(&mut self, src: &mut Source) -> &mut Self;
    fn time_propagate(&mut self, secs: f64, src: &mut Source) -> &mut Self;
}

/// Wrapper for CEO source
///
/// # Examples
///
/// ```
/// use gicsdom::ceo::Source;
/// let mut src = Source::new(1,25.5,401);
/// src.build("V",vec![0.0],vec![0.0],vec![0.0]);
/// ```
pub struct Source {
    pub _c_: source,
    /// The number of sources
    size: i32,
    /// The diameter of the entrance pupil [m] 
    pub pupil_size: f64,
    /// The sampling of the entrance pupil [px]
    pub pupil_sampling: i32,
    pub _wfe_rms: Vec<f32>,
}
impl Source {
    /// Creates and empty `Source`
    pub fn empty() -> Source {
        Source {
            _c_: unsafe { mem::zeroed() },
            size: 0,
            pupil_size: 0.0,
            pupil_sampling: 0,
            _wfe_rms: vec![],
        }
    }
    /// Creates a new `Source` with the arguments:
    ///
    /// * `pupil_size` - the diameter of the entrance pupil [m] 
    /// * `pupil_sampling` - the sampling of the entrance pupil [px] 
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
    /// Sets the `Source` parameters:
    ///
    /// * `band` - the photometric band: Vs, V, R, I, J, H, K or R+I
    /// * `zenith` - the zenith angle [rd]
    /// * `azimuth` - the azimuth angle [rd]
    /// * `magnitude` - the magnitude at the specified photometric band
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
    /// Returns the `Source` wavelength [m]
    pub fn wavelength(&mut self) -> f64 {
        unsafe { self._c_.wavelength() as f64 }
    }
    /// Sets the `Source` full width at half maximum in un-binned detector pixel
    pub fn set_fwhm(&mut self, value: f64) {
        self._c_.fwhm = value as f32;
    }
    /// Set the pupil rotation angle [degree]
    pub fn rotate_rays(&mut self, angle: f64) {
        self._c_.rays.rot_angle = angle;
    }
    /// Copies the optical path difference from ray tracing into the wavefront phase argument, this usually takes place after ray tracing to the exit pupil
    pub fn xpupil(&mut self) -> &mut Self {
        unsafe {
            self._c_.wavefront.reset();
            self._c_.opd2phase();
        }
        self
    }
    /// Returns the wavefront error root mean square [m]
    pub fn wfe_rms(&mut self) -> Vec<f32> {
        unsafe {
            self._c_.wavefront.rms(self._wfe_rms.as_mut_ptr());
        }
        self._wfe_rms.clone()
    }
    /// Returns the wavefront error root mean square [m]x10^-`exp`
    pub fn wfe_rms_10e(&mut self, exp: i32) -> Vec<f32> {
        unsafe {
            self._c_.wavefront.rms(self._wfe_rms.as_mut_ptr());
        }
        self._wfe_rms
            .iter()
            .map(|x| x * 10_f32.powi(-exp))
            .collect()
    }
    /// Returns the x and y gradient of the wavefront in average over each of the GMT segments
    pub fn segments_gradients(&mut self) -> Vec<Vec<f32>> {
        let mut sxy: Vec<Vec<f32>> = vec![vec![0.; 7 * self.size as usize]; 2];
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
    /// Returns the x and y gradient of the wavefront in average over each lenslet of a `n_lenslet`x`n_lenslet` array, the gradients are saved in `Centroiding`
    pub fn lenslet_gradients(&mut self, n_lenslet: i32, _lenslet_size: f64, data: &mut Centroiding) {
        let lenslet_size = self.pupil_size/n_lenslet as f64;
        unsafe {
            if data.n_valid_lenslet < data.n_lenslet_total {
                self._c_.wavefront.finite_difference1(
                    data.__mut_ceo__().0.d__cx,
                    data.__mut_ceo__().0.d__cy,
                    n_lenslet,
                    lenslet_size as f32,
                    data.__mut_ceo__().1,
                );
            } else {
                self._c_.wavefront.finite_difference(
                    data.__mut_ceo__().0.d__cx,
                    data.__mut_ceo__().0.d__cy,
                    n_lenslet,
                    lenslet_size as f32,
                );
            }
        }
    }
    /// Resets the rays and the wavefront to their original state
    pub fn reset(&mut self) {
        unsafe {
            self._c_.wavefront.reset();
            self._c_.reset_rays();
        }
    }
    /// Updates the `zenith` and `azimuth` of the `Source`
    pub fn update(&mut self, mut zenith: Vec<f64>, mut azimuth: Vec<f64>) {
        unsafe {
            self._c_.update_directions(
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                zenith.len() as i32,
            );
        }
    }
    /// Adds `phase` to the `Source` wavefront
    pub fn add(&mut self, phase: &mut CuFloat) -> &mut Self {
        unsafe {
            self._c_.wavefront.add_phase(1.0, phase.as_mut_ptr());
        }
        self
    }
    /// Propagates a `Source` through a `system` that implements the `Propagation` trait
    pub fn through<T: Propagation>(&mut self, system: &mut T) -> &mut Self {
        system.propagate(self);
        self
    }
}
impl Drop for Source {
    /// Frees CEO memory before dropping `Source`
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
