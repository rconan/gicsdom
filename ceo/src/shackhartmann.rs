//!
//! # CEO shackhartmann wrapper
//!
//! Provides a structure `ShackHartmann` that is a wrapper for [CEO](https://github.com/rconan/CEO) shackhartmann C++ structure.
//! `ShackHartmann<M: Model>` is instantiated and initialized with the `SHACKHARTMANN<M: Model>` builder where `Model` is either type `Geometric` of `Diffractive`
//!
//! # Examples
//!
//! ```
//! use ceo::ceo;
//! // Creates a gmt instance with default parameters
//! let mut wfs = ceo!(SHACKHARTMANN:Geometric);
//! ```

use super::ceo_bindings::{geometricShackHartmann, shackHartmann};
use super::{cu::Single, Builder, Cu, Mask, Propagation, Source, SOURCE};
use std::{f32, mem};

pub type Geometric = geometricShackHartmann;
pub type Diffractive = shackHartmann;

pub trait Model: Clone {
    fn build(
        &mut self,
        n_side_lenslet: i32,
        d: f32,
        n_sensor: i32,
        n_px_lenslet: i32,
        osf: i32,
        n_px: i32,
        b: i32,
    );
    fn get_c_as_mut_ptr(&mut self) -> *mut f32;
    fn drop(&mut self);
}
impl Model for Geometric {
    fn build(
        &mut self,
        n_side_lenslet: i32,
        d: f32,
        n_sensor: i32,
        _n_px_lenslet: i32,
        _osf: i32,
        _n_px: i32,
        _b: i32,
    ) {
        unsafe {
            self.setup(n_side_lenslet, d, n_sensor);
        }
    }
    fn get_c_as_mut_ptr(&mut self) -> *mut f32 {
        self.data_proc.d__c
    }
    fn drop(&mut self) {
        unsafe { self.cleanup() };
    }
}
impl Model for Diffractive {
    fn drop(&mut self) {
        unsafe {
            self.cleanup();
        }
    }
    fn get_c_as_mut_ptr(&mut self) -> *mut f32 {
        self.data_proc.d__c
    }
    fn build(
        &mut self,
        n_side_lenslet: i32,
        d: f32,
        n_sensor: i32,
        n_px_lenslet: i32,
        osf: i32,
        n_px: i32,
        b: i32,
    ) {
        unsafe {
            self.setup(n_side_lenslet, n_px_lenslet, d, osf, n_px, b, n_sensor);
        }
    }
}

// n_side_lenslet, n_px_lenslet, d
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct LensletArray(pub usize, pub usize, pub f64);
impl Default for LensletArray {
    fn default() -> Self {
        LensletArray(1, 511, 25.5)
    }
}
// n_px_framelet, n_px_imagelet, osf
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct Detector(pub usize, pub Option<usize>, pub Option<usize>);
impl Default for Detector {
    fn default() -> Self {
        Detector(512, None, None)
    }
}
/// `ShackHartmann` builder
///
/// Default properties:
///  - n_sensor: 1
///  - lenslet_array:
///    - n_lenslet: 1
///    - n_px_lenslet: 511px
///    - lenslet_pitch: 25.5m
///  - detector:
///    - n_px_framelet: 512px
///    - n_px_imagelet: None[512px]
///    - osf: None[2]
///
/// # Examples
///
/// ```
/// use ceo::{Builder,  SHACKHARTMANN, Geometric};
/// let mut wfs = SHACKHARTMANN::<Geometric>::new().build();
/// ```
#[derive(Debug, Clone)]
pub struct SHACKHARTMANN<T: Model> {
    pub n_sensor: usize,
    pub lenslet_array: LensletArray,
    pub detector: Detector,
    marker: std::marker::PhantomData<T>,
}
impl<T: Model> Default for SHACKHARTMANN<T> {
    fn default() -> Self {
        SHACKHARTMANN {
            n_sensor: 1,
            lenslet_array: LensletArray::default(),
            detector: Detector::default(),
            marker: std::marker::PhantomData,
        }
    }
}
impl<T: Model> SHACKHARTMANN<T> {
    pub fn set_n_sensor(self, n_sensor: usize) -> Self {
        Self { n_sensor, ..self }
    }
    pub fn set_lenslet_array(self, n_side_lenslet: usize, n_px_lenslet: usize, d: f64) -> Self {
        Self {
            lenslet_array: LensletArray(n_side_lenslet, n_px_lenslet, d),
            ..self
        }
    }
    pub fn guide_stars(&self) -> SOURCE {
        let LensletArray(n_side_lenslet, n_px_lenslet, d) = self.lenslet_array;
        SOURCE::new()
            .set_size(self.n_sensor)
            .set_pupil_size(d * n_side_lenslet as f64)
            .set_pupil_sampling(n_px_lenslet * n_side_lenslet + 1)
    }
}
impl<T: Model> Builder for SHACKHARTMANN<T> {
    type Component = ShackHartmann<T>;
    fn build(self) -> ShackHartmann<T> {
        let LensletArray(n_side_lenslet, n_px_lenslet, d) = self.lenslet_array;
        let mut wfs = ShackHartmann::<T> {
            _c_: unsafe { std::mem::zeroed() },
            n_side_lenslet: n_side_lenslet as i32,
            n_px_lenslet: n_px_lenslet as i32,
            d,
            n_sensor: self.n_sensor as i32,
            n_centroids: 0,
            centroids: super::Cu::vector(
                (n_side_lenslet * n_side_lenslet * 2 * self.n_sensor) as usize,
            ),
        };
        let Detector(n_px_framelet, n_px_imagelet, osf) = self.detector;
        let n_px = match n_px_imagelet {
            Some(n_px_imagelet) => n_px_imagelet,
            None => n_px_framelet,
        };
        let b = n_px / n_px_framelet;
        let o = match osf {
            Some(osf) => osf,
            None => 2,
        };
        wfs.n_centroids = wfs.n_side_lenslet * wfs.n_side_lenslet * 2 * wfs.n_sensor;
        wfs._c_.build(
            wfs.n_side_lenslet,
            wfs.d as f32,
            wfs.n_sensor,
            wfs.n_px_lenslet,
            o as i32,
            n_px as i32,
            b as i32,
        );
        wfs.centroids.from_ptr(wfs._c_.get_c_as_mut_ptr());
        wfs
    }
}
/// `ShackHartmann` "SH48" builder for GMT AGWS model
///
/// Default properties:
///  - n_sensor: 4
///  - lenslet_array:
///    - n_lenslet: 48
///    - n_px_lenslet: 16px
///    - lenslet_pitch: 25.5m/48
///  - detector:
///    - n_px_framelet: 8px
///    - n_px_imagelet: Some(24px)
///    - osf: Some(2)
///
/// # Examples
///
/// ```
/// use ceo::{Builder, SH48, Geometric};
/// let mut wfs = SH48::<Geometric>::new().build();
/// ```
#[derive(Debug, Clone)]
pub struct SH48<T: Model> {
    pub n_sensor: usize,
    pub lenslet_array: LensletArray,
    pub detector: Detector,
    marker: std::marker::PhantomData<T>,
}
impl<T: Model> Default for SH48<T> {
    fn default() -> Self {
        SH48 {
            n_sensor: 4,
            lenslet_array: LensletArray(48, 16, 25.5 / 48.0),
            detector: Detector(8, Some(24), Some(2)),
            marker: std::marker::PhantomData,
        }
    }
}
impl<T: Model> SH48<T> {
    pub fn set_n_sensor(self, n_sensor: usize) -> Self {
        Self { n_sensor, ..self }
    }
    pub fn guide_stars(&self) -> SOURCE {
        let LensletArray(n_side_lenslet, n_px_lenslet, d) = self.lenslet_array;
        SOURCE::new()
            .set_size(self.n_sensor)
            .set_pupil_size(d * n_side_lenslet as f64)
            .set_pupil_sampling(n_px_lenslet * n_side_lenslet + 1)
    }
}
impl<T: Model> Builder for SH48<T> {
    type Component = ShackHartmann<T>;
    fn build(self) -> ShackHartmann<T> {
        let LensletArray(n_side_lenslet, n_px_lenslet, d) = self.lenslet_array;
        let mut wfs = ShackHartmann::<T> {
            _c_: unsafe { std::mem::zeroed() },
            n_side_lenslet: n_side_lenslet as i32,
            n_px_lenslet: n_px_lenslet as i32,
            d,
            n_sensor: self.n_sensor as i32,
            n_centroids: 0,
            centroids: super::Cu::vector(
                (n_side_lenslet * n_side_lenslet * 2 * self.n_sensor) as usize,
            ),
        };
        let Detector(n_px_framelet, n_px_imagelet, osf) = self.detector;
        let n_px = match n_px_imagelet {
            Some(n_px_imagelet) => n_px_imagelet,
            None => n_px_framelet,
        };
        let b = n_px / n_px_framelet;
        let o = match osf {
            Some(osf) => osf,
            None => 2,
        };
        wfs.n_centroids = wfs.n_side_lenslet * wfs.n_side_lenslet * 2 * wfs.n_sensor;
        wfs._c_.build(
            wfs.n_side_lenslet,
            wfs.d as f32,
            wfs.n_sensor,
            wfs.n_px_lenslet,
            o as i32,
            n_px as i32,
            b as i32,
        );
        wfs.centroids.from_ptr(wfs._c_.get_c_as_mut_ptr());
        wfs
    }
}
/// shackhartmann wrapper
pub struct ShackHartmann<S: Model> {
    pub _c_: S,
    /// The size of the square lenslet array
    pub n_side_lenslet: i32,
    /// The number of pixel per lenslet in the telescope pupil
    pub n_px_lenslet: i32,
    /// The lenslet array pitch [m]
    pub d: f64,
    /// The number of WFS
    pub n_sensor: i32,
    /// The total number of centroids
    pub n_centroids: i32,
    /// The centroids
    pub centroids: Cu<Single>,
}
impl<S: Model> ShackHartmann<S> {
    /// Creates a new `ShackHartmann` as either `Geometric` or `Diffractive` type
    ///
    /// * `n_sensor` - the number of WFS
    /// * `n_side_lenslet` - the size of the square lenslet array
    /// * `n_px_lenslet` - the number of pixel per lenslet in the telescope pupil
    /// * `d` - the lenslet pitch [m]
    pub fn new(n_sensor: i32, n_side_lenslet: i32, n_px_lenslet: i32, d: f64) -> ShackHartmann<S> {
        ShackHartmann {
            _c_: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
            n_centroids: 0,
            centroids: Cu::vector((n_side_lenslet * n_side_lenslet * 2 * n_sensor) as usize),
        }
    }
    pub fn guide_stars(&self) -> SOURCE {
        SOURCE::new()
            .set_size(self.n_sensor as usize)
            .set_pupil_size(self.d * self.n_side_lenslet as f64)
            .set_pupil_sampling((self.n_px_lenslet * self.n_side_lenslet + 1) as usize)
            .set_band("R")
    }
    pub fn drop(&mut self) {
        self._c_.drop();
    }
}
impl ShackHartmann<Geometric> {
    /// Initializes the `ShackHartmann` WFS
    pub fn build(&mut self) -> &mut Self {
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        unsafe {
            self._c_
                .setup(self.n_side_lenslet, self.d as f32, self.n_sensor);
            self.centroids.from_ptr(self._c_.data_proc.d__c);
        }
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
    /// Calibrates the `ShackHartmann` WFS reference slopes and valid lenslets
    pub fn calibrate(&mut self, src: &mut Source, threshold: f64) -> &mut Self {
        unsafe {
            self._c_.calibrate(src.as_raw_mut_ptr(), threshold as f32);
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        unsafe {
            self._c_.process();
        }
        self
    }
    pub fn get_data(&mut self) -> Cu<Single> {
        let m = self._c_.valid_lenslet.nnz as usize * 2usize;
        let mut data: Cu<Single> = Cu::vector(m);
        data.malloc();
        unsafe {
            self._c_.get_valid_slopes(data.as_ptr());
        }
        data
    }
    pub fn filter(&mut self, lenslet_mask: &mut Mask) -> Cu<Single> {
        let m = lenslet_mask.nnz() as usize * 2usize;
        let mut data: Cu<Single> = Cu::vector(m);
        data.malloc();
        unsafe {
            self._c_
                .masked_slopes(data.as_ptr(), lenslet_mask.as_mut_prt());
        }
        data
    }
    pub fn fold_into(&mut self, data: &mut Cu<Single>, lenslet_mask: &mut Mask) {
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
    pub fn lenset_mask(&mut self) -> Cu<Single> {
        let mut mask: Cu<Single> =
            Cu::vector((self.n_side_lenslet * self.n_side_lenslet * self.n_sensor) as usize);
        mask.from_ptr(self._c_.valid_lenslet.f);
        mask
    }
    pub fn lenlet_flux(&mut self) -> Cu<Single> {
        let mut flux: Cu<Single> =
            Cu::vector((self.n_side_lenslet * self.n_side_lenslet * self.n_sensor) as usize);
        flux.from_ptr(self._c_.data_proc.d__mass);
        flux
    }
    pub fn as_raw_mut_ptr(&mut self) -> &mut Geometric {
        &mut self._c_
    }
}
impl<S: Model> Drop for ShackHartmann<S> {
    fn drop(&mut self) {
        self.drop();
    }
}
impl Propagation for ShackHartmann<Geometric> {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.propagate(src.as_raw_mut_ptr());
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
            self._c_.setup(
                self.n_side_lenslet,
                self.n_px_lenslet,
                self.d as f32,
                o,
                n_px,
                b,
                self.n_sensor,
            );
            self.centroids.from_ptr(self._c_.data_proc.d__c);
        }
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
            self._c_.calibrate(src.as_raw_mut_ptr(), threshold as f32);
            self._c_.camera.reset();
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        unsafe {
            self._c_.process();
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
impl Propagation for ShackHartmann<Diffractive> {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.propagate(src.as_raw_mut_ptr());
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}

impl From<ShackHartmann<Geometric>> for Source {
    fn from(item: ShackHartmann<Geometric>) -> Self {
        item.new_guide_stars()
    }
}
impl From<ShackHartmann<Diffractive>> for Source {
    fn from(item: ShackHartmann<Diffractive>) -> Self {
        item.new_guide_stars()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shack_hartmann_geometric_new() {
        use crate::GMT;
        let mut wfs = SHACKHARTMANN::<Geometric>::new()
            .set_n_sensor(1)
            .set_lenslet_array(48, 16, 25.5 / 48f64)
            .build();
        let mut src = SOURCE::new().set_pupil_sampling(48 * 16 + 1).build();
        let mut gmt = GMT::new().build();
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }

    #[test]
    fn shack_hartmann_geometric_new_with_macro() {
        let mut wfs = crate::ceo!(
            SHACKHARTMANN:Geometric,
            set_n_sensor = [1],
            set_lenslet_array = [48, 16, 25.5 / 48f64]
        );
        let mut src = crate::ceo!(SOURCE, set_pupil_sampling = [48 * 16 + 1]);
        let mut gmt = crate::ceo!(GMT);
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }

    /*
    #[test]
    fn shack_hartmann_diffractive_new() {
        use crate::Builder;
        use element::*;
        let mut wfs = CEO::<SHACKHARTMANN<Diffractive>>::new()
            .set_n_sensor(1)
            .set_lenslet_array(48, 16, 25.5 / 48f64)
            .set_detector(8, Some(24), None)
            .build();
        let mut src = CEO::<SOURCE>::new().set_pupil_sampling(48 * 16 + 1).build();
        let mut gmt = CEO::<GMT>::new().build();
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }

    #[test]
    fn shack_hartmann_diffractive_new_with_macro() {
        use element::*;
        let mut wfs = crate::ceo!(
            SHACKHARTMANN<Diffractive>,
            set_n_sensor = [1],
            set_lenslet_array = [48, 16, 25.5 / 48f64],
            set_detector = [8, Some(24), None]
        );
        let mut src = crate::ceo!(SOURCE, set_pupil_sampling = [48 * 16 + 1]);
        let mut gmt = crate::ceo!(GMT);
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }
    */
}
