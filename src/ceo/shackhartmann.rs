use std::ffi::CString;
use std::{f32, mem};

use super::ceo_bindings::{dev2host,geometricShackHartmann,shackHartmann};
use super::Propagation;
use super::Source;

pub struct Geometric_ShackHartmann {
    _c_: geometricShackHartmann,
    pub n_side_lenslet: i32,
    pub n_px_lenslet: i32,
    pub d: f64,
    pub n_sensor: i32,
    pub n_centroids: i32,
    pub centroids: Vec<f32>,
}
pub struct Diffractive_ShackHartmann {
    _c_: shackHartmann,
    n_side_lenslet: i32,
    n_px_lenslet: i32,
    d: f64,
    n_sensor: i32,
    pub n_centroids: i32,
    pub centroids: Vec<f32>,
}
impl Geometric_ShackHartmann {
    pub fn new(
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        d: f64,
    ) -> Geometric_ShackHartmann {
        Geometric_ShackHartmann {
            _c_: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
            n_centroids: 0,
            centroids: Vec::new(),
        }
    }
    pub fn build(&mut self) -> &mut Self {
        unsafe {
            self._c_
                .setup(self.n_side_lenslet, self.d as f32, self.n_sensor);
        }
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        self.centroids = vec![0.0; self.n_centroids as usize]; //Vec::with_capacity(self.n_centroids as usize);
        self
    }
    pub fn guide_star_args(&self) -> (i32, f64, i32) {
        (
            self.n_sensor,
            self.d * self.n_side_lenslet as f64,
            self.n_px_lenslet * self.n_side_lenslet + 1,
        )
    }
    pub fn new_guide_stars(&self) -> Source {
        Source::new(
            self.n_sensor,
            self.d * self.n_side_lenslet as f64,
            self.n_px_lenslet * self.n_side_lenslet + 1,
        )
    }
    pub fn calibrate(&mut self, src: &mut Source, threshold: f64) -> &mut Self {
        unsafe {
            self._c_.calibrate(&mut src._c_, threshold as f32);
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        self.centroids = vec![0.0; self.n_centroids as usize];
        //let mut c = vec![0.0;self.n_centroids as usize];
        unsafe {
            self._c_.process();
            dev2host(
                self.centroids.as_mut_ptr(),
                self._c_.data_proc.d__c,
                self.n_centroids as i32,
            );
            self._c_.reset();
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.reset();
        }
        self
    }
}
impl Drop for Geometric_ShackHartmann {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Geometric_ShackHartmann {
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
impl Diffractive_ShackHartmann {
    pub fn new(
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        d: f64,
    ) -> Diffractive_ShackHartmann {
        Diffractive_ShackHartmann {
            _c_: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
            n_centroids: 0,
            centroids: Vec::new(),
        }
    }
    pub fn build(
        &mut self,
        n_px_framelet: i32,
        n_px_imagelet: Option<i32>,
        osf: Option<i32>,
    ) -> &mut Self {
        let n_px = match n_px_imagelet {
            Some(n_px_imagelet) => n_px_imagelet,
            None => n_px_framelet,
        };
        let b = n_px / n_px_framelet;
        let o = match osf {
            Some(osf) => osf,
            None => 2,
        };
        unsafe {
            self._c_.setup(
                self.n_side_lenslet,
                self.n_px_lenslet,
                self.d as f32,
                o,
                n_px,
                b,
                self.n_sensor,
            );
        }
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        self.centroids = vec![0.0; self.n_centroids as usize]; //Vec::with_capacity(self.n_centroids as usize);
        self
    }
    pub fn new_guide_stars(&self) -> Source {
        Source::new(
            self.n_sensor,
            self.d * self.n_side_lenslet as f64,
            self.n_px_lenslet * self.n_side_lenslet + 1,
        )
    }
    pub fn calibrate(&mut self, src: &mut Source, threshold: f64) -> &mut Self {
        unsafe {
            self._c_.calibrate(&mut src._c_, threshold as f32);
            self._c_.camera.reset();
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        self.centroids = vec![0.0; self.n_centroids as usize];
        //let mut c = vec![0.0;self.n_centroids as usize];
        unsafe {
            self._c_.process();
            dev2host(
                self.centroids.as_mut_ptr(),
                self._c_.data_proc.d__c,
                self.n_centroids as i32,
            );
            self._c_.camera.reset();
        }
        self
    }
    pub fn readout(
        &mut self,
        exposure_time: f32,
        readout_noise_rms: f32,
        n_background_photon: f32,
        noise_factor: f32,
    ) -> &mut Self {
        unsafe {
            self._c_.camera.readout1(
                exposure_time,
                readout_noise_rms,
                n_background_photon,
                noise_factor,
            );
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.camera.reset();
        }
        self
    }
}
impl Drop for Diffractive_ShackHartmann {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Diffractive_ShackHartmann {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.propagate(&mut src._c_);
        }
        self
    }
    fn time_propagate(&mut self, secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}
pub type GeometricShackHartmann = Geometric_ShackHartmann;
pub type ShackHartmann = Diffractive_ShackHartmann;
