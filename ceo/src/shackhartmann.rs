use std::{f32, mem};

use super::ceo_bindings::{geometricShackHartmann, shackHartmann};
use super::{Cu, Mask, Propagation, Source};

pub type Geometric = geometricShackHartmann;
pub type Diffractive = shackHartmann;

enum Model {
    None,
    Geometric,
    Diffractive,
}

pub struct ShackHartmann<S> {
    _c_: std::marker::PhantomData<S>,
    _c_geometric: Geometric,
    _c_diffractive: Diffractive,
    pub n_side_lenslet: i32,
    pub n_px_lenslet: i32,
    pub d: f64,
    pub n_sensor: i32,
    pub n_centroids: i32,
    pub centroids: Cu<f32>,
    model: Model,
}
impl<S> ShackHartmann<S> {
    pub fn new(n_sensor: i32, n_side_lenslet: i32, n_px_lenslet: i32, d: f64) -> ShackHartmann<S> {
        ShackHartmann {
            _c_: std::marker::PhantomData,
            _c_geometric: unsafe { mem::zeroed() },
            _c_diffractive: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
            n_centroids: 0,
            centroids: Cu::vector((n_side_lenslet * n_side_lenslet * 2 * n_sensor) as usize),
            model: Model::None,
        }
    }
}
impl ShackHartmann<Geometric> {
    pub fn build(&mut self) -> &mut Self {
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        unsafe {
            self._c_geometric
                .setup(self.n_side_lenslet, self.d as f32, self.n_sensor);
            self.centroids.from_ptr(self._c_geometric.data_proc.d__c);
        }
        self.model = Model::Geometric;
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
            self._c_geometric.calibrate(&mut src._c_, threshold as f32);
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        unsafe {
            self._c_geometric.process();
        }
        self
    }
    pub fn get_data(&mut self) -> Cu<f32> {
        let m = self._c_geometric.valid_lenslet.nnz as usize * 2usize;
        let mut data: Cu<f32> = Cu::vector(m);
        data.malloc();
        unsafe {
            self._c_geometric.get_valid_slopes(data.as_ptr());
        }
        data
    }
    pub fn filter(&mut self, lenslet_mask: &mut Mask) -> Cu<f32> {
        let m = lenslet_mask.nnz() as usize * 2usize;
        let mut data: Cu<f32> = Cu::vector(m);
        data.malloc();
        unsafe {
            self._c_geometric
                .masked_slopes(data.as_ptr(), lenslet_mask.as_mut_prt());
        }
        data
    }
    pub fn fold_into(&mut self, data: &mut Cu<f32>, lenslet_mask: &mut Mask) {
        unsafe {
            self._c_geometric
                .folded_slopes(data.as_ptr(), lenslet_mask.as_mut_prt());
        }
    }
    pub fn n_valid_lenslet(&mut self) -> usize {
        self._c_geometric.valid_lenslet.nnz as usize
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_geometric.reset();
        }
        self
    }
    pub fn lenset_mask(&mut self) -> Cu<f32> {
        let mut mask: Cu<f32> =
            Cu::vector((self.n_side_lenslet * self.n_side_lenslet * self.n_sensor) as usize);
        mask.from_ptr(self._c_geometric.valid_lenslet.f);
        mask
    }
    pub fn lenlet_flux(&mut self) -> Cu<f32> {
        let mut flux: Cu<f32> =
            Cu::vector((self.n_side_lenslet * self.n_side_lenslet * self.n_sensor) as usize);
        flux.from_ptr(self._c_geometric.data_proc.d__mass);
        flux
    }
}
impl<S> Drop for ShackHartmann<S> {
    fn drop(&mut self) {
        match self.model {
            Model::None => (),
            Model::Geometric => {
                unsafe { self._c_geometric.cleanup() };
            }
            Model::Diffractive => {
                unsafe { self._c_diffractive.cleanup() };
            }
        }
    }
}
impl Propagation for ShackHartmann<Geometric> {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_geometric.propagate(&mut src._c_);
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}
impl ShackHartmann<Diffractive> {
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
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        unsafe {
            self._c_diffractive.setup(
                self.n_side_lenslet,
                self.n_px_lenslet,
                self.d as f32,
                o,
                n_px,
                b,
                self.n_sensor,
            );
            self.centroids.from_ptr(self._c_diffractive.data_proc.d__c);
        }
        self.model = Model::Diffractive;
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
            self._c_diffractive
                .calibrate(&mut src._c_, threshold as f32);
            self._c_diffractive.camera.reset();
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        unsafe {
            self._c_diffractive.process();
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
            self._c_diffractive.camera.readout1(
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
            self._c_diffractive.camera.reset();
        }
        self
    }
}
impl Propagation for ShackHartmann<Diffractive> {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_diffractive.propagate(&mut src._c_);
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}
