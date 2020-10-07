use std::{f32, mem};

use super::ceo_bindings::{geometricShackHartmann, shackHartmann};
use super::{
    element, CeoElement, CeoError, Cu, Detector, LensletArray, Mask, Propagation, Source, CEO,
};

pub type Geometric = geometricShackHartmann;
pub type Diffractive = shackHartmann;

pub trait Model {
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

/// Shack-Hartmann wavefront sensor model
pub struct ShackHartmann<S: Model> {
    _c_: S,
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
    pub centroids: Cu<f32>,
}
impl CEO<element::SHACKHARTMANN> {
    pub fn new() -> CEO<element::SHACKHARTMANN> {
        CEO {
            element: CeoElement::Shackhartmann {
                n_sensor: 1,
                lenslet_array: LensletArray(1, 511, 25.5),
                detector: Detector(512, None, None),
            },
            t: std::marker::PhantomData,
        }
    }
    pub fn set_n_sensor(mut self, n_sensor: usize) -> Self {
        if let CeoElement::Shackhartmann {
            n_sensor: _,
            lenslet_array,
            detector,
        } = self.element
        {
            self.element = CeoElement::Shackhartmann {
                n_sensor: n_sensor,
                lenslet_array,
                detector,
            }
        }
        self
    }
    pub fn set_lenslet_array(mut self, n_side_lenslet: usize, n_px_lenslet: usize, d: f64) -> Self {
        if let CeoElement::Shackhartmann {
            n_sensor,
            lenslet_array: _,
            detector,
        } = self.element
        {
            self.element = CeoElement::Shackhartmann {
                n_sensor,
                lenslet_array: LensletArray(n_side_lenslet, n_px_lenslet, d),
                detector,
            }
        }
        self
    }
    pub fn set_detector(
        mut self,
        n_px_framelet: usize,
        n_px_imagelet: Option<usize>,
        osf: Option<usize>,
    ) -> Self {
        if let CeoElement::Shackhartmann {
            n_sensor,
            lenslet_array,
            detector: _,
        } = self.element
        {
            self.element = CeoElement::Shackhartmann {
                n_sensor: n_sensor,
                lenslet_array,
                detector: Detector(n_px_framelet, n_px_imagelet, osf),
            }
        }
        self
    }
    pub fn build<T: Model>(self) -> Result<ShackHartmann<T>, CeoError<element::SHACKHARTMANN>> {
        match self.element {
            CeoElement::Shackhartmann {
                n_sensor,
                lenslet_array,
                detector,
            } => {
                let LensletArray(n_side_lenslet, n_px_lenslet, d) = lenslet_array;
                let mut wfs = ShackHartmann::<T> {
                    _c_: unsafe { mem::zeroed() },
                    n_side_lenslet: n_side_lenslet as i32,
                    n_px_lenslet: n_px_lenslet as i32,
                    d,
                    n_sensor: n_sensor as i32,
                    n_centroids: 0,
                    centroids: Cu::vector(
                        (n_side_lenslet * n_side_lenslet * 2 * n_sensor) as usize,
                    ),
                };
                let Detector(n_px_framelet, n_px_imagelet, osf) = detector;
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
                Ok(wfs)
            }
            _ => Err(CeoError(element::SHACKHARTMANN)),
        }
    }
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
            self._c_.calibrate(&mut src._c_, threshold as f32);
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        unsafe {
            self._c_.process();
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
impl<S: Model> Drop for ShackHartmann<S> {
    fn drop(&mut self) {
        self.drop();
    }
}
impl Propagation for ShackHartmann<Geometric> {
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
            self._c_.calibrate(&mut src._c_, threshold as f32);
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
            self._c_.propagate(&mut src._c_);
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
        use element::*;
        let mut wfs = CEO::<SHACKHARTMANN>::new()
            .set_n_sensor(1)
            .set_lenslet_array(48, 16, 25.5 / 48f64)
            .build::<Geometric>()
            .unwrap();
        let mut src = CEO::<SOURCE>::new()
            .set_pupil_sampling(48 * 16 + 1)
            .build()
            .unwrap();
        let mut gmt = CEO::<GMT>::new().build().unwrap();
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }

    #[test]
    fn shack_hartmann_geometric_new_with_macro() {
        use element::*;
        let mut wfs = crate::ceo!(
            SHACKHARTMANN: Geometric,
            set_n_sensor = [1],
            set_lenslet_array = [48, 16, 25.5 / 48f64]
        );
        let mut src = crate::ceo!(SOURCE, set_pupil_sampling = [48 * 16 + 1]);
        let mut gmt = crate::ceo!(GMT);
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }

    #[test]
    fn shack_hartmann_diffractive_new() {
        use element::*;
        let mut wfs = CEO::<SHACKHARTMANN>::new()
            .set_n_sensor(1)
            .set_lenslet_array(48, 16, 25.5 / 48f64)
            .set_detector(8, Some(24), None)
            .build::<Diffractive>()
            .unwrap();
        let mut src = CEO::<SOURCE>::new()
            .set_pupil_sampling(48 * 16 + 1)
            .build()
            .unwrap();
        let mut gmt = CEO::<GMT>::new().build().unwrap();
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }

    #[test]
    fn shack_hartmann_diffractive_new_with_macro() {
        use element::*;
        let mut wfs = crate::ceo!(
            SHACKHARTMANN: Diffractive,
            set_n_sensor = [1],
            set_lenslet_array = [48, 16, 25.5 / 48f64],
            set_detector = [8, Some(24), None]
        );
        let mut src = crate::ceo!(SOURCE, set_pupil_sampling = [48 * 16 + 1]);
        let mut gmt = crate::ceo!(GMT);
        src.through(&mut gmt).xpupil().through(&mut wfs);
        println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    }
}
