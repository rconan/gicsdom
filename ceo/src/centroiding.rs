use std::mem;

use super::ceo_bindings::{centroiding, dev2host, host2dev_char, mask};
use super::Imaging;

/// Wrapper for CEO centroiding
pub struct Centroiding {
    _c_: centroiding,
    _c_mask_: mask,
    /// The total number of lenslets
    pub n_lenslet_total: u32,
    /// The number of centroids
    pub n_centroids: u32,
    /// The centroid units, default: 1 (pixel)
    pub units: f32,
    flux: Vec<f32>,
    /// The valid lenslet mask
    pub valid_lenslets: Vec<i8>,
    /// The number of valid lenslet
    pub n_valid_lenslet: u32,
    /// The centroids
    pub centroids: Vec<f32>,
}

impl Centroiding {
    /// Creates a new `Centroiding`
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
    /// Sets the `Centroiding` parameters:
    ///
    /// * `n_lenslet` - the size of the square lenslet array
    /// * `data_units` - the centroids units
    pub fn build(&mut self, n_lenslet: u32, data_units: Option<f64>) -> &mut Self {
        self.n_lenslet_total = n_lenslet * n_lenslet;
        self.n_centroids = 2 * self.n_lenslet_total;
        self.n_valid_lenslet = self.n_lenslet_total;
        unsafe {
            self._c_.setup(n_lenslet as i32, 1);
            self._c_mask_.setup(self.n_lenslet_total as i32);
        }
        self.flux = vec![0.0; self.n_lenslet_total as usize];
        self.centroids = vec![0.0; self.n_centroids as usize];
        self.units = data_units.or(Some(1f64)).unwrap() as f32;
        self
    }
    /// Computes the `centroids` from the `sensor` image; optionally, a `Centroiding` `reference` may be provided that offsets the `centroids` and sets the `valid_lenslets`
    pub fn process(&mut self, sensor: &Imaging, reference: Option<&Centroiding>) -> &mut Self {
        if reference.is_none() {
            unsafe {
                self._c_.get_data2(
                    sensor.__ceo__().d__frame,
                    sensor.__ceo__().N_PX_CAMERA,
                    0.0,
                    0.0,
                    self.units,
                );
            }
        } else {
            let r = reference.unwrap();
            assert_eq!(self.n_lenslet_total, r.n_lenslet_total);
            unsafe {
                self._c_.get_data3(
                    sensor.__ceo__().d__frame,
                    sensor.__ceo__().N_PX_CAMERA,
                    r._c_.d__cx,
                    r._c_.d__cy,
                    self.units,
                    r._c_mask_.m,
                );
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
        }
        self
    }
    /// grabs the `centroids` from the GPU
    pub fn grab(&mut self) -> &mut Self {
        unsafe {
            dev2host(
                self.centroids.as_mut_ptr(),
                self._c_.d__c,
                self.n_centroids as i32,
            );
        }
        self
    }
    /// returns the valid `centroids` i.e. the `centroids` that corresponds to a non-zero entry in the `valid_lenslets` mask; if `some_valid_lenslet` is given, then it supersedes any preset `valid_lenset`
    pub fn valids(&self, some_valid_lenslets: Option<&Vec<i8>>) -> Vec<f32> {
        let valid_lenslets = some_valid_lenslets.or(Some(&self.valid_lenslets)).unwrap();
        assert_eq!(self.n_lenslet_total, valid_lenslets.len() as u32);
        let n = valid_lenslets.iter().fold(0u32, |a, x| a + (*x as u32)) as usize;
        let mut valid_centroids: Vec<f32> = vec![0f32; 2 * n];
        let mut l = 0;
        for (k, v) in valid_lenslets.iter().enumerate() {
            if *v > 0 {
                valid_centroids[l] = self.centroids[k];
                valid_centroids[l + n] = self.centroids[k + valid_lenslets.len()];
                l += 1;
            }
        }
        return valid_centroids;
    }
    /// returns the flux of each lenslet
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
    /// returns the sum of the flux of all the lenslets
    pub fn integrated_flux(&mut self) -> f32 {
        self.lenslet_flux().iter().sum()
    }
    /// Computes the valid lenslets and return the number of valid lenslets; the valid lenslets are computed based on the maximum flux threshold or a given valid lenslets mask
    pub fn set_valid_lenslets(
        &mut self,
        some_flux_threshold: Option<f64>,
        some_valid_lenslets: Option<Vec<i8>>,
    ) -> u32 {
        if some_flux_threshold.is_some() {
            let lenslet_flux = self.lenslet_flux();
            let lenslet_flux_max = lenslet_flux.iter().cloned().fold(0.0, f32::max);
            let threshold_flux = lenslet_flux_max * some_flux_threshold.unwrap() as f32;
            self.valid_lenslets = lenslet_flux
                .iter()
                .map(|x| if x >= &threshold_flux { 1i8 } else { 0i8 })
                .collect();
        }
        if some_valid_lenslets.is_some() {
            self.valid_lenslets = some_valid_lenslets.unwrap();
        }
        unsafe {
            host2dev_char(
                self._c_.lenslet_mask,
                self.valid_lenslets.as_mut_ptr(),
                self.n_lenslet_total as i32,
            );
            self._c_mask_.reset();
            self._c_mask_
                .add1(self._c_.lenslet_mask, self.n_lenslet_total as i32);
        }
        self.n_valid_lenslet = self
            .valid_lenslets
            .iter()
            .fold(0u32, |a, x| a + (*x as u32));
        self.n_valid_lenslet
    }
    pub fn __ceo__(&mut self) -> (&centroiding, &mask) {
        (&self._c_, &self._c_mask_)
    }
    pub fn __mut_ceo__(&mut self) -> (&mut centroiding, &mut mask) {
        (&mut self._c_, &mut self._c_mask_)
    }
}
impl Drop for Centroiding {
    /// Frees CEO memory before dropping `Centroiding`
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
            self._c_mask_.cleanup();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Gmt,Source,Conversion};

    #[test]
    fn centroiding_sim() {
        let pupil_size = 25.5f64;
        let n_side_lenslet = 48;
        let n_px_lenslet = 16;
        let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;
        let lenslet_size = (pupil_size / n_side_lenslet as f64) as f32;
        let mut gmt = Gmt::new();
        gmt.build(0, None);
        let mut src = Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0f32], vec![0f32], vec![18f32]);
        src.set_fwhm(3f64);
        let mut sensor = Imaging::new();
        sensor.build(1, n_side_lenslet, n_px_lenslet, 2, 24, 3);
        let p = sensor.pixel_scale(&mut src) as f64;

        let mut cog0 = Centroiding::new();
        cog0.build(n_side_lenslet as u32, None);
        src.through(&mut gmt).xpupil().through(&mut sensor);
        let nv = cog0
            .process(&sensor, None)
            .set_valid_lenslets(Some(0.5), None);
        println!("Valid lenslet #: {}",nv);

        let mut cog = Centroiding::new();
        cog.build(n_side_lenslet as u32, Some(p))
            .set_valid_lenslets(None, Some(cog0.valid_lenslets.clone()));
        src.through(&mut gmt).xpupil().lenslet_gradients(
            n_side_lenslet,
            lenslet_size as f64,
            &mut cog,
        );
        let s0 = cog.grab().valids(None);

        let m2_rbm = vec![vec![0f64, 0f64, 0f64, 1e-6, 1e-6, 0f64]; 7];
        gmt.update(None, Some(&m2_rbm), None, None);
        sensor.reset();
        src.through(&mut gmt).xpupil().through(&mut sensor);
        let c = cog.process(&sensor, Some(&cog0)).grab().valids(None);
        src.lenslet_gradients(n_side_lenslet, lenslet_size as f64, &mut cog);
        let s = cog
            .grab()
            .valids(None)
            .iter()
            .zip(s0.iter())
            .map(|x| x.0 - x.1)
            .collect::<Vec<f32>>();

        let e = ((c
            .iter()
            .zip(s.iter())
            .map(|x| (x.0 - x.1).powi(2))
            .sum::<f32>()
            / nv as f32)
            .sqrt() as f64)
            .to_mas();
        println!("Centroid error: {}mas", e);
        assert!(e<5f64);
    }
}
