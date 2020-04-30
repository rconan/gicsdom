use std::mem;

use super::ceo_bindings::{dev2host, imaging};
use super::Propagation;
use super::Source;

#[derive(Copy, Clone)]
/// A square lenslet array
pub struct LensletArray {
    /// The number of lenslet per side
    pub n_side_lenslet: i32,
    /// Dimension [m] of one lenslet
    pub lenslet_size: f64,
}
#[derive(Copy, Clone)]
/// Detector noise specifications
pub struct NoiseDataSheet {
    /// Read-out noise rms
    pub rms_read_out_noise: f64,
    /// Number of background photons
    pub n_background_photon: f64,
    /// Noise factor
    pub noise_factor: f64,
}

/// An optical imager with a detector
///
/// The optical imager is a square lenslet array which focal plane lies on the detector.
/// The detector continuously integrates the images formed on the detector until it is explicitly reset
pub struct Imaging {
    _c_: imaging,
    dft_osf: i32,
}
impl Imaging {
    /// Creates a new `Imaging`
    pub fn new() -> Imaging {
        Imaging {
            _c_: unsafe { mem::zeroed() },
            dft_osf: 1,
        }
    }
    /// Set `Imaging` parameters
    ///
    /// * `n_sensor` - the number of `Imaging` sensor
    /// * `n_side_lenslet` - the size of the square lenslet array
    /// * `n_px_lenslet` - the number of pixel per lenslet, for a total resolution of (`n_side_lenslet`X`n_px_lenslet`+1)^2
    /// * `dft_osf` - the discrete Fourier transform oversampling factor
    /// * `n_px_imagelet` - the sampling of a lenslet focal plane image
    /// * `binning` - binning factor of a imagelet
    pub fn build(
        &mut self,
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        dft_osf: i32,
        n_px_imagelet: i32,
        binning: i32,
    ) -> &mut Self {
        unsafe {
            self._c_.setup3(
                n_px_lenslet,
                n_side_lenslet,
                dft_osf,
                n_px_imagelet,
                binning,
                n_sensor,
            );
        }
        self.dft_osf = dft_osf;
        self
    }
    pub fn __ceo__(&self) -> &imaging {
        &self._c_
    }
    /// Returns the frame from the GPU
    pub fn frame_transfer(&mut self, frame: &mut Vec<f32>) -> &mut Self {
        unsafe {
            dev2host(
                frame.as_mut_ptr(),
                self._c_.d__frame,
                self.resolution() * self.resolution() * self._c_.N_SOURCE,
            );
        }
        self
    }
    /// Resets the detector frame to zero
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.reset();
        }
        self
    }
    /// Returns the detector resolution
    pub fn resolution(&self) -> i32 {
        self._c_.N_PX_CAMERA * self._c_.N_SIDE_LENSLET
    }
    /// Return the number of frames that have been summed since the last reset
    pub fn n_frame(&self) -> u32 {
        self._c_.N_FRAME as u32
    }
    /// Reads out the detector by adding noise to the frame if a `NoiseDataSheet` is passed and the intensity is scaled according to the detector `exposure` time [s]
    pub fn readout(
        &mut self,
        exposure: f64,
        detector_noise_properties: Option<NoiseDataSheet>,
    ) -> &mut Self {
        detector_noise_properties.map_or_else(
            || (),
            |p| unsafe {
                self._c_.readout1(
                    exposure as f32,
                    p.rms_read_out_noise as f32,
                    p.n_background_photon as f32,
                    p.noise_factor as f32,
                );
            },
        );
        self
    }
    /// Sets the pixel scale
    pub fn set_pixel_scale(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.pixel_scale = (src.wavelength() / src.pupil_size) as f32
                * (self._c_.N_SIDE_LENSLET * self._c_.BIN_IMAGE / self.dft_osf) as f32;
        }
        self
    }
    /// Sets the detector pointing direction
    pub fn set_pointing(&mut self, mut zen: Vec<f32>, mut azi: Vec<f32>) -> &mut Self {
        unsafe {
            self._c_.absolute_pointing = 1;
            self._c_
                .set_pointing_direction(zen.as_mut_ptr(), azi.as_mut_ptr());
        }
        self
    }
}
impl Drop for Imaging {
    /// Frees CEO memory before dropping `Imaging`
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Imaging {
    /// Fourier propagates the wavefront to the focal plane onto the detector
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.propagate(&mut src._c_);
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}
