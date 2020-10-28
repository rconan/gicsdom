//!
//! # CEO wrapper crate
//!
//!The CEO wrapper is the interface to [CEO CUDA API](https://github.com/rconan/CEO).
//! The simplest method to build CEO element is to use the [`ceo!`][macro] macro builder
//!
//! [macro]: macro.ceo.html

use std::{error::Error, f32, f64, fmt, mem};

pub mod atmosphere;
pub mod calibrations;
pub mod centroiding;
pub mod ceo_bindings;
pub mod cu;
pub mod fwhm;
pub mod gmt;
pub mod imaging;
pub mod lmmse;
pub mod pssn;
pub mod shackhartmann;
pub mod source;

#[doc(inline)]
pub use self::atmosphere::Atmosphere;
#[doc(inline)]
pub use self::calibrations::Calibration;
#[doc(inline)]
pub use self::centroiding::Centroiding;
#[doc(inline)]
pub use self::cu::Cu;
#[doc(inline)]
pub use self::fwhm::Fwhm;
#[doc(inline)]
pub use self::gmt::Gmt;
#[doc(inline)]
pub use self::imaging::Imaging;
#[doc(inline)]
pub use self::lmmse::LinearMinimumMeanSquareError;
#[doc(inline)]
pub use self::pssn::PSSn;
#[doc(inline)]
pub use self::shackhartmann::ShackHartmann;
#[doc(inline)]
pub use self::source::Propagation;
#[doc(inline)]
pub use self::source::Source;
#[doc(hidden)]
pub use ceo_bindings::{geqrf, gpu_double, gpu_float, mask, ormqr, set_device};

pub type GeometricShackHartmann = ShackHartmann<shackhartmann::Geometric>;

/// CEO macro builder
///
/// One macro to rule them all, one macro to find them, one macro to bring them all and in the darkness bind them all
///
/// # Examples
///
///  * GMT
///
/// ```
/// use ceo::{ceo, element::*};
/// let gmt = ceo!(GMT, set_m1_n_mode = [27], set_m2_n_mode = [123]);
/// ```
///
///  * Geometric Shack-Hartmann
///
/// ```
/// use ceo::{ceo, element::*, shackhartmann::Geometric};
/// let mut wfs = ceo!(
///     SHACKHARTMANN: Geometric,
///     set_n_sensor = [1],
///     set_lenslet_array = [48, 16, 25.5 / 48f64]
/// );
/// let mut src = ceo!(SOURCE, set_pupil_sampling = [48 * 16 + 1]);
/// let mut gmt = ceo!(GMT);
/// src.through(&mut gmt).xpupil().through(&mut wfs);
/// println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
/// ```
///
///  * Diffractive Shack-Hartmann
///
/// ```
/// use ceo::{ceo, element::*, shackhartmann::Diffractive};
/// let mut wfs = ceo!(
///     SHACKHARTMANN: Diffractive,
///     set_n_sensor = [1],
///     set_lenslet_array = [48, 16, 25.5 / 48f64],
///     set_detector = [8, Some(24), None]
/// );
/// let mut src = ceo!(SOURCE, set_pupil_sampling = [48 * 16 + 1]);
/// let mut gmt = ceo!(GMT);
/// src.through(&mut gmt).xpupil().through(&mut wfs);
/// println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
/// ```
#[macro_export]
macro_rules! ceo {
    ($element:ty) => {
        CEO::<$element>::new().build()
    };
    ($element:ty:$model:ty) => {
        CEO::<$element>::new().build::<$model>()
    };
    ($element:ty, $($arg:ident = [$($val:expr),+]),*) => {
        CEO::<$element>::new()$(.$arg($($val),+))*.build()
    };
    ($element:ty:$model:ty, $($arg:ident = [$($val:expr),+]),*) => {
        CEO::<$element>::new()$(.$arg($($val),+))*.build::<$model>()
    };
}
/*
macro_rules! gmt {
    ($($arg:ident = $val:expr),*) => {
        ceo!(element::GMT,$($arg = $val),*)
    };
}
*/
#[derive(Debug)]
pub struct CeoError<T>(T);
impl<T: std::fmt::Debug> Error for CeoError<T> {}
impl<T: std::fmt::Debug> fmt::Display for CeoError<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "CEO {:?} builder has failed!", self.0)
    }
}

/// CEO builder type trait
///
/// Only structures in the [`element`][element] module implement the trait
///
/// [element]: element/index.html
pub trait CEOType: Clone + Default {}
/// CEO builder pattern
///
/// `CEO` is a generic builder pattern for all CEO elements.
/// It will accept only the structures of the [`element`][element] module that implements the [`CEOtype`][ceotype] trait.
///
/// [element]: element/index.html
/// [ceotype]: trait.CEOType.html
#[derive(Debug, Clone)]
pub struct CEO<T: CEOType> {
    args: T,
}
pub trait CEOInit<T: CEOType> {
    fn new() -> CEO<T>;
}
impl<T: CEOType> CEOInit<T> for CEO<T> {
    fn new() -> Self {
        CEO {
            args: Default::default(),
        }
    }
}
pub trait CEOWFS {
    fn build(&self) -> ShackHartmann<shackhartmann::Geometric>;
    fn get_n_data(&self) -> usize;
}
pub mod element {
    use super::CEOType;
    #[doc(hidden)]
    #[derive(Debug, Clone)]
    pub struct Mirror {
        pub mode_type: String,
        pub n_mode: usize,
    }
    /// [`CEO`](../struct.CEO.html#impl) [`Gmt`](../struct.Gmt.html) builder type
    #[derive(Debug, Clone)]
    pub struct GMT {
        pub m1: Mirror,
        pub m2: Mirror,
    }
    /// Default properties:
    ///  * M1:
    ///    * mode type : "bending modes"
    ///    * \# mode    : 0
    ///  * M2:
    ///    * mode type : "Karhunen-Loeve"
    ///    * \# mode    : 0
    impl Default for GMT {
        fn default() -> Self {
            GMT {
                m1: Mirror {
                    mode_type: "bending modes".into(),
                    n_mode: 0,
                },
                m2: Mirror {
                    mode_type: "Karhunen-Loeve".into(),
                    n_mode: 0,
                },
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    #[derive(Debug, Clone)]
    /// [`CEO`](../struct.CEO.html#impl-4) [`Source`](../struct.Source.html) builder type
    pub struct SOURCE {
        pub size: usize,
        pub pupil_size: f64,
        pub pupil_sampling: usize,
        pub band: String,
        pub zenith: Vec<f32>,
        pub azimuth: Vec<f32>,
        pub magnitude: Vec<f32>,
    }
    /// Default properties:
    ///  * size             : 1
    ///  * pupil size       : 25.5m
    ///  * pupil sampling   : 512px
    ///  * photometric band : Vs (500nm)
    ///  * zenith           : 0degree
    ///  * azimuth          : 0degree
    ///  * magnitude        : 0
    impl Default for SOURCE {
        fn default() -> Self {
            SOURCE {
                size: 1,
                pupil_size: 25.5,
                pupil_sampling: 512,
                band: "Vs".into(),
                zenith: vec![0f32],
                azimuth: vec![0f32],
                magnitude: vec![0f32],
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    #[derive(Debug, Clone)]
    /// [`CEO`](../struct.CEO.html#impl-3) specialized [`Source`](../struct.Source.html) builder type
    pub struct FIELDDELAUNAY21 {
        pub size: usize,
        pub pupil_size: f64,
        pub pupil_sampling: usize,
        pub band: String,
        pub zenith: Vec<f32>,
        pub azimuth: Vec<f32>,
        pub magnitude: Vec<f32>,
    }
    use serde::{Deserialize, Serialize};
    #[derive(Debug, Deserialize, Serialize, Default)]
    struct GlaoField {
        pub zenith_arcmin: Vec<f32>,
        pub azimuth_degree: Vec<f32>,
    }
    impl Default for FIELDDELAUNAY21 {
        fn default() -> Self {
            use super::Conversion;
            use serde_pickle as pickle;
            use std::fs::File;
            let field_reader = File::open("ceo/fielddelaunay21.pkl").expect("File not found!");
            let field: GlaoField = pickle::from_reader(field_reader).expect("File loading failed!");
            let n_src = field.zenith_arcmin.len();
            FIELDDELAUNAY21 {
                size: n_src,
                pupil_size: 25.5,
                pupil_sampling: 512,
                band: "Vs".into(),
                zenith: field
                    .zenith_arcmin
                    .iter()
                    .map(|x| x.from_arcmin())
                    .collect::<Vec<f32>>(),
                azimuth: field
                    .azimuth_degree
                    .iter()
                    .map(|x| x.to_radians())
                    .collect::<Vec<f32>>(),
                magnitude: vec![0f32; n_src],
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
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
    pub trait ShWfs {
        fn build<T: super::shackhartmann::Model>(
            &self,
            n_sensor: usize,
            lenslet_array: LensletArray,
            detector: Detector,
        ) -> super::ShackHartmann<T> {
            let LensletArray(n_side_lenslet, n_px_lenslet, d) = lenslet_array;
            let mut wfs = super::ShackHartmann::<T> {
                _c_: unsafe { std::mem::zeroed() },
                n_side_lenslet: n_side_lenslet as i32,
                n_px_lenslet: n_px_lenslet as i32,
                d,
                n_sensor: n_sensor as i32,
                n_centroids: 0,
                centroids: super::Cu::vector(
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
            wfs
        }
        fn guide_stars(&self, n_sensor: usize, lenslet_array: LensletArray) -> super::CEO<SOURCE> {
            let LensletArray(n_side_lenslet, n_px_lenslet, d) = lenslet_array;
            super::CEO::<SOURCE>::new()
                .set_size(n_sensor)
                .set_pupil_size(d * n_side_lenslet as f64)
                .set_pupil_sampling(n_px_lenslet * n_side_lenslet + 1)
        }
    }
    #[derive(Debug, Clone)]
    /// [`CEO`](../struct.CEO.html#impl-2) [`ShackHartmann`](../struct.ShackHartmann.html) builder type
    pub struct SHACKHARTMANN {
        pub n_sensor: usize,
        pub lenslet_array: LensletArray,
        pub detector: Detector,
    }
    impl Default for SHACKHARTMANN {
        fn default() -> Self {
            SHACKHARTMANN {
                n_sensor: 1,
                lenslet_array: LensletArray::default(),
                detector: Detector::default(),
            }
        }
    }
    impl ShWfs for SHACKHARTMANN {}
    #[derive(Debug, Clone)]
    /// [`CEO`](../struct.CEO.html#impl-7) specialized [`ShackHartmann`](../struct.ShackHartmann.html) builder type
    pub struct SH48 {
        pub n_sensor: usize,
        pub lenslet_array: LensletArray,
        pub detector: Detector,
    }
    impl Default for SH48 {
        fn default() -> Self {
            SH48 {
                n_sensor: 4,
                lenslet_array: LensletArray(48, 16, 25.5 / 48.0),
                detector: Detector(8, Some(24), Some(2)),
            }
        }
    }
    impl ShWfs for SH48 {}
    // ---------------------------------------------------------------------------------------------
    #[derive(Debug, Clone)]
    /// [`CEO`](../struct.CEO.html#impl-1) [`PSSn`](../struct.PSSn.html) builder type
    pub struct PSSN {
        pub r0_at_zenith: f64,
        pub oscale: f64,
        pub zenith_angle: f64,
    }
    /// Default properties:
    ///  * r0           : 16cm
    ///  * L0           : 25m
    ///  * zenith angle : 30 degrees
    impl Default for PSSN {
        fn default() -> Self {
            PSSN {
                r0_at_zenith: 0.16,
                oscale: 25.0,
                zenith_angle: 30_f64.to_radians(),
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    #[derive(Debug, Clone)]
    #[doc(hidden)]
    pub struct TurbulenceProfile {
        pub n_layer: usize,
        pub altitude: Vec<f32>,
        pub xi0: Vec<f32>,
        pub wind_speed: Vec<f32>,
        pub wind_direction: Vec<f32>,
    }
    impl Default for TurbulenceProfile {
        fn default() -> Self {
            TurbulenceProfile {
                n_layer: 7,
                altitude: [25.0, 275.0, 425.0, 1_250.0, 4_000.0, 8_000.0, 13_000.0].to_vec(),
                xi0: [0.1257, 0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751].to_vec(),
                wind_speed: [5.6540, 5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187].to_vec(),
                wind_direction: [0.0136, 0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462].to_vec(),
            }
        }
    }
    #[derive(Debug, Clone)]
    #[doc(hidden)]
    pub struct RayTracing {
        pub width: f32,
        pub n_width_px: i32,
        pub field_size: f32,
        pub duration: f32,
        pub filepath: Option<String>,
        pub n_duration: Option<i32>,
    }
    /// [`CEO`](../struct.CEO.html#impl-6) [`Atmosphere`](../struct.Atmosphere.html) builder type
    #[derive(Debug, Clone)]
    pub struct ATMOSPHERE {
        pub r0_at_zenith: f64,
        pub oscale: f64,
        pub zenith_angle: f64,
        pub turbulence: TurbulenceProfile,
        pub ray_tracing: Option<RayTracing>,
    }
    /// Default properties:
    ///  * r0           : 16cm
    ///  * L0           : 25m
    ///  * zenith angle : 30 degrees
    ///  * turbulence profile:
    ///    * n_layer        : 7
    ///    * altitude       : [25.0, 275.0, 425.0, 1250.0, 4000.0, 8000.0, 13000.0] m
    ///    * xi0            : [0.1257, 0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751]
    ///    * wind speed     : [5.6540, 5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187] m/s
    ///    * wind direction : [0.0136, 0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462] rd
    /// * ray tracing : none
    impl Default for ATMOSPHERE {
        fn default() -> Self {
            ATMOSPHERE {
                r0_at_zenith: 0.16,
                oscale: 25.,
                zenith_angle: 30_f64.to_radians(),
                turbulence: TurbulenceProfile::default(),
                ray_tracing: None,
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    #[derive(Debug, Clone)]
    pub struct LMMSE {
        pub atm: super::CEO<ATMOSPHERE>,
        pub guide_star: super::CEO<SOURCE>,
        pub mmse_star: super::CEO<SOURCE>,
        pub fov_diameter: Option<f64>,
        pub n_side_lenslet: usize,
        pub solver_id: String,
        pub wavefront_osf: usize,
    }
    impl Default for LMMSE {
        fn default() -> Self {
            LMMSE {
                atm: super::CEO::<ATMOSPHERE>::new(),
                guide_star: super::CEO::<SOURCE>::new(),
                mmse_star: super::CEO::<SOURCE>::new(),
                fov_diameter: None,
                n_side_lenslet: 0,
                solver_id: "MINRES".to_owned(),
                wavefront_osf: 1,
            }
        }
    }
    // ---------------------------------------------------------------------------------------------
    macro_rules! impl_ceotype {
        ($($element:ty),+) => {
            $(impl CEOType for $element {})+
        };
    }
    impl_ceotype!(
        GMT,
        SOURCE,
        SHACKHARTMANN,
        PSSN,
        FIELDDELAUNAY21,
        ATMOSPHERE,
        SH48,
        LMMSE
    );
}

pub trait Conversion<T> {
    fn from_arcmin(self) -> T;
    fn from_arcsec(self) -> T;
    fn from_mas(self) -> T;
    fn to_arcmin(self) -> T;
    fn to_arcsec(self) -> T;
    fn to_mas(self) -> T;
}
macro_rules! impl_conversion {
    ($($name:ty),+) => {
        $(impl Conversion<$name> for $name {
            /// Converts angle in arcminute to radian
            fn from_arcmin(self) -> $name {
                self.to_radians() / 60.
            }
            /// Converts angle in arcsecond to radian
            fn from_arcsec(self) -> $name {
                self.from_arcmin() / 60.
            }
            /// Converts angle in milli-arcsecond to radian
            fn from_mas(self) -> $name {
                self.from_arcsec() * 1e-3
            }
            /// Converts angle in radian to arcminute
            fn to_arcmin(self) -> $name {
                60.0 * self.to_degrees()
            }
            /// Converts angle in radian to arcsecond
            fn to_arcsec(self) -> $name {
                60.0 * self.to_arcmin()
            }
            /// Converts angle in radian to mill-arcsecond
            fn to_mas(self) -> $name {
                1e3 * self.to_arcsec()
            }
        })+
    };
}
impl_conversion!(f64, f32);

pub fn set_gpu(id: i32) {
    unsafe {
        set_device(id);
    }
}

use cu::Single;
#[derive(Clone, Debug)]
pub struct Mask {
    _c_: mask,
}
impl Mask {
    pub fn new() -> Self {
        Mask {
            _c_: unsafe { mem::zeroed() },
        }
    }
    pub fn build(&mut self, n_el: usize) -> &mut Self {
        unsafe { self._c_.setup(n_el as i32) }
        self
    }
    pub fn filter(&mut self, f: &mut Cu<Single>) -> &mut Self {
        unsafe {
            self._c_.alter(f.as_ptr());
            self._c_.set_index();
        }
        self
    }
    pub fn nnz(&self) -> usize {
        self._c_.nnz as usize
    }
    pub fn as_mut_prt(&mut self) -> *mut mask {
        &mut self._c_
    }
    pub fn as_raw_mut_ptr(&mut self) -> &mut mask {
        &mut self._c_
    }
}
impl Default for Mask {
    fn default() -> Self {
        Self::new()
    }
}

pub struct CuFloat {
    _c_: gpu_float,
    host_data: Vec<f32>,
}
impl CuFloat {
    pub fn new() -> CuFloat {
        CuFloat {
            _c_: unsafe { mem::zeroed() },
            host_data: vec![],
        }
    }
    pub fn malloc(&mut self, len: usize) -> &mut Self {
        unsafe {
            self._c_.setup1(len as i32);
            self._c_.dev_malloc();
        }
        self
    }
    pub fn into(&mut self, host_data: &mut Vec<f32>) -> &mut Self {
        unsafe {
            self._c_.host_data = host_data.as_mut_ptr();
            self._c_.host2dev();
        }
        self
    }
    pub fn from(&mut self) -> Vec<f32> {
        self.host_data = vec![0f32; self._c_.N as usize];
        unsafe {
            self._c_.host_data = self.host_data.as_mut_ptr();
            self._c_.dev2host()
        }
        self.host_data.clone()
    }
    pub fn as_mut_ptr(&mut self) -> *mut f32 {
        self._c_.dev_data
    }
    pub fn mv(&mut self, y: &mut CuFloat, x: &mut CuFloat) -> &mut Self {
        unsafe {
            self._c_.mv(&mut y._c_, &mut x._c_);
        }
        self
    }
    pub fn qr(&mut self, m: i32) -> &mut Self {
        unsafe {
            self._c_.qr(m);
        }
        self
    }
    pub fn qr_solve(&mut self, x: &mut CuFloat, b: &mut CuFloat) -> &mut Self {
        unsafe {
            self._c_.qr_solve(&mut x._c_, &mut b._c_);
        }
        self
    }
}
impl Drop for CuFloat {
    fn drop(&mut self) {
        unsafe { self._c_.free_dev() }
    }
}

pub fn qr(tau: &mut CuFloat, a: &mut CuFloat, m: i32, n: i32) {
    unsafe { geqrf(tau._c_.dev_data, a._c_.dev_data, m, n) }
}
pub fn qtb(b: &mut CuFloat, m: i32, a: &mut CuFloat, tau: &mut CuFloat, n: i32) {
    unsafe { ormqr(b._c_.dev_data, m, a._c_.dev_data, tau._c_.dev_data, n) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn cufloat_mv() {
        let mut a = CuFloat::new();
        a.malloc(9)
            .into(&mut (0..9).map(|x| x as f32).collect::<Vec<f32>>());
        let mut x = CuFloat::new();
        x.malloc(3).into(&mut vec![1f32; 3]);
        let mut y = CuFloat::new();
        y.malloc(3);
        a.mv(&mut y, &mut x);
        println!("R: {:?}", y.from());
        assert_eq!(y.from(), vec![9f32, 12f32, 15f32]);
    }

    #[test]
    fn cufloat_qr() {
        let mut a = CuFloat::new();
        a.malloc(12)
            .into(&mut vec![
                1f32, 4f32, 2f32, 1f32, 2f32, 5f32, 1f32, 1f32, 3f32, 6f32, 1f32, 1f32,
            ])
            .qr(4);
        println!("a: {:?}", a.from());
        let mut b = CuFloat::new();
        b.malloc(4).into(&mut vec![6f32, 15f32, 4f32, 3f32]);
        let mut x = CuFloat::new();
        x.malloc(3);
        a.qr_solve(&mut x, &mut b);
        println!("x: {:?}", x.from());
    }

    #[test]
    fn cufloat_bigqr() {
        let m = 5000;
        let n = 700;
        let p = m * n;
        let mut rng = rand::thread_rng();
        let mut a = CuFloat::new();
        a.malloc(p).into(
            &mut (0..p)
                .map(|_| rng.gen_range(-100, 100) as f32)
                .collect::<Vec<f32>>(),
        );
        //println!("a: {:?}",a.from());
        let mut b = CuFloat::new();
        b.malloc(m);
        {
            let mut x = CuFloat::new();
            x.malloc(n).into(&mut vec![1f32; n]);
            a.mv(&mut b, &mut x);
        }
        //println!("b: {:?}",b.from());
        a.qr(m as i32);
        //println!("a: {:?}",a.from());
        let mut x = CuFloat::new();
        x.malloc(n);
        a.qr_solve(&mut x, &mut b);
        let sx = x.from().iter().sum::<f32>();
        println!("sum of x: {:?}", sx);
        assert!((sx - (n as f32)).abs() < 1e-6);
    }
}
