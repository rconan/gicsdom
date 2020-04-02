use std::mem;

use super::ceo_bindings::{centroiding, dev2host, host2dev_char, mask};
use super::Imaging;

pub struct Centroiding {
    _c_: centroiding,
    _c_mask_: mask,
    pub n_lenslet_total: u32,
    pub n_centroids: u32,
    pub units: f32,
    flux: Vec<f32>,
    pub valid_lenslets: Vec<i8>,
    pub n_valid_lenslet: u32,
    pub centroids: Vec<f32>,
}

impl Centroiding {
    pub fn new() -> Centroiding {
        Centroiding {
            _c_: unsafe { mem::zeroed() },
            _c_mask_: unsafe { mem::zeroed() },
            n_lenslet_total: 0u32,
            n_centroids: 0u32,
            units: 1f32,
            flux: vec![],
            valid_lenslets: vec![],
            n_valid_lenslet: 0u32,
            centroids: vec![],
        }
    }
    pub fn build(&mut self, n_lenslet: u32, n_source: u32) -> &mut Self {
        self.n_lenslet_total = n_lenslet * n_lenslet * n_source;
        self.n_centroids = 2 * self.n_lenslet_total;
        self.n_valid_lenslet = self.n_lenslet_total;
        unsafe {
            self._c_.setup(n_lenslet as i32, n_source as i32);
            self._c_mask_.setup(self.n_lenslet_total as i32);
        }
        self.flux = vec![0.0; self.n_lenslet_total as usize];
        self.centroids = vec![0.0; self.n_centroids as usize];
        self
    }
    pub fn process(
        &mut self,
        sensor: &Imaging,
        reference: Option<&Centroiding>,
    ) -> Option<Vec<f32>> {
        if reference.is_none() {
            unsafe {
                self._c_
                    .get_data1(sensor.__ceo__().d__frame, sensor.__ceo__().N_PX_CAMERA);
            }
            return None;
        } else {
            let r = reference.unwrap();
            unsafe {
                self._c_.get_data3(
                    sensor.__ceo__().d__frame,
                    sensor.__ceo__().N_PX_CAMERA,
                    r._c_.d__cx,
                    r._c_.d__cy,
                    self.units,
                    r._c_mask_.m,
                );
                dev2host(
                    self.centroids.as_mut_ptr(),
                    self._c_.d__c,
                    self.n_centroids as i32,
                );
            }
            let mut valid_centroids: Vec<f32> = vec![];
            if r.n_valid_lenslet < self.n_lenslet_total {
                let n = r.n_valid_lenslet as usize;
                valid_centroids = vec![0f32; 2 * n];
                let mut l = 0;
                for (k, v) in r.valid_lenslets.iter().enumerate() {
                    if *v > 0 {
                        valid_centroids[l] = self.centroids[k];
                        valid_centroids[l + n] = self.centroids[k + r.n_lenslet_total as usize];
                        l += 1;
                    }
                }
            }
            /*
            let n = r.n_lenslet_total as usize;
            let cx = &self.centroids[..n];
            let cy = &self.centroids[n..];
            let mut  vcx: Vec<f32> = r
                .valid_lenslets
                .iter()
                .zip(cx.iter())
                .filter(|x| x.0.is_positive())
                .map(|x| *x.1)
                .collect();
            let mut vcy: Vec<f32> = r
                .valid_lenslets
                .iter()
                .zip(cy.iter())
                .filter(|x| x.0.is_positive())
                .map(|x| *x.1)
                .collect();
            vcx.append(&mut vcy);
            */
            return Some(valid_centroids.clone());
        }
    }
    pub fn lenslet_flux(&mut self) -> &Vec<f32> {
        unsafe {
            dev2host(
                self.flux.as_mut_ptr(),
                self._c_.d__mass,
                self.n_lenslet_total as i32,
            );
        }
        &self.flux
    }
    pub fn integrated_flux(&mut self) -> f32 {
        self.lenslet_flux().iter().sum()
    }
    pub fn set_valid_lenslets(&mut self, flux_threshold: f64) -> u32 {
        let lenslet_flux = self.lenslet_flux();
        let lenslet_flux_max = lenslet_flux.iter().cloned().fold(0.0, f32::max);
        let threshold_flux = lenslet_flux_max * flux_threshold as f32;
        self.valid_lenslets = lenslet_flux
            .iter()
            .map(|x| if x > &threshold_flux { 1i8 } else { 0i8 })
            .collect();
        unsafe {
            host2dev_char(
                self._c_.lenslet_mask,
                self.valid_lenslets.as_mut_ptr(),
                self.n_lenslet_total as i32,
            );
        }
        self.n_valid_lenslet = self
            .valid_lenslets
            .iter()
            .fold(0u32, |a, x| a + (*x as u32));
        self.n_valid_lenslet
    }
}
impl Drop for Centroiding {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
            self._c_mask_.cleanup();
        }
    }
}
