use std::{f32, mem};

use super::ceo_bindings::{dev2host, geometricShackHartmann, shackHartmann};
use super::{Cu, Mask, Propagation, Source};

pub struct GeometricShackHartmann {
    _c_: geometricShackHartmann,
    pub n_side_lenslet: i32,
    pub n_px_lenslet: i32,
    pub d: f64,
    pub n_sensor: i32,
    pub n_centroids: i32,
    pub centroids: Cu<f32>,
}
pub struct DiffractiveShackHartmann {
    _c_: shackHartmann,
    n_side_lenslet: i32,
    n_px_lenslet: i32,
    d: f64,
    n_sensor: i32,
    pub n_centroids: i32,
    pub centroids: Vec<f32>,
}
impl GeometricShackHartmann {
    pub fn new(
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        d: f64,
    ) -> GeometricShackHartmann {
        GeometricShackHartmann {
            _c_: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
            n_centroids: 0,
            centroids: Cu::vector((n_side_lenslet * n_side_lenslet * 2 * n_sensor) as usize),
        }
    }
    pub fn build(&mut self) -> &mut Self {
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        unsafe {
            self._c_
                .setup(self.n_side_lenslet, self.d as f32, self.n_sensor);
            self.centroids.from_ptr(self._c_.data_proc.d__c);
        }
        //self.centroids = vec![0.0; self.n_centroids as usize]; //Vec::with_capacity(self.n_centroids as usize);
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
        //self.centroids = vec![0.0; self.n_centroids as usize];
        //let mut c = vec![0.0;self.n_centroids as usize];
        unsafe {
            self._c_.process();
            /*
            dev2host(
                self.centroids.as_mut_ptr(),
                self._c_.data_proc.d__c,
                self.n_centroids as i32,
            );
            */
            //self._c_.reset();
        }
        self
    }
    pub fn get_data(&mut self) -> Cu<f32> {
        let m = self._c_.valid_lenslet.nnz as usize * 2usize;
        let mut data: Cu<f32> = Cu::vector(m);
        data.malloc();
        unsafe {
            self._c_.get_valid_slopes(data.as_ptr());
        }
        data
    }
    pub fn filter(&mut self, lenslet_mask: &mut Mask) -> Cu<f32> {
        let m = lenslet_mask.nnz() as usize * 2usize;
        let mut data: Cu<f32> = Cu::vector(m);
        data.malloc();
        unsafe {
            self._c_
                .masked_slopes(data.as_ptr(), lenslet_mask.as_mut_prt());
        }
        data
    }
    pub fn fold_into(&mut self, data: &mut Cu<f32>, lenslet_mask: &mut Mask) {
        unsafe {
            self._c_
                .folded_slopes(data.as_ptr(), lenslet_mask.as_mut_prt());
        }
    }
    pub fn n_valid_lenslet(&mut self) -> usize {
        self._c_.valid_lenslet.nnz as usize
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.reset();
        }
        self
    }
    pub fn lenset_mask(&mut self) -> Cu<f32> {
        let mut mask: Cu<f32> =
            Cu::vector((self.n_side_lenslet * self.n_side_lenslet * self.n_sensor) as usize);
        mask.from_ptr(self._c_.valid_lenslet.f);
        mask
    }
    pub fn lenlet_flux(&mut self) -> Cu<f32> {
        let mut flux: Cu<f32> =
            Cu::vector((self.n_side_lenslet * self.n_side_lenslet * self.n_sensor) as usize);
        flux.from_ptr(self._c_.data_proc.d__mass);
        flux
    }
}
impl Drop for GeometricShackHartmann {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for GeometricShackHartmann {
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
impl DiffractiveShackHartmann {
    pub fn new(
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        d: f64,
    ) -> DiffractiveShackHartmann {
        DiffractiveShackHartmann {
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
impl Drop for DiffractiveShackHartmann {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for DiffractiveShackHartmann {
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
pub type ShackHartmann = DiffractiveShackHartmann;
