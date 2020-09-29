use std::{f32, f64, mem};

pub mod atmosphere;
//pub mod calibrations;
pub mod centroiding;
pub mod ceo_bindings;
pub mod cu;
pub mod fwhm;
pub mod gmt;
pub mod imaging;
pub mod pssn;
pub mod shackhartmann;
pub mod source;

pub use self::atmosphere::Atmosphere;
//pub use self::calibrations::Calibration;
pub use self::centroiding::Centroiding;
pub use self::cu::Cu;
pub use self::fwhm::Fwhm;
pub use self::gmt::Gmt;
pub use self::imaging::{Imaging, LensletArray};
pub use self::pssn::PSSn;
pub use self::shackhartmann::ShackHartmann;
pub use self::source::Propagation;
pub use self::source::Source;
pub use ceo_bindings::{geqrf, gpu_double, gpu_float, mask, ormqr, set_device};

pub type GeometricShackHartmann = ShackHartmann<shackhartmann::Geometric>;

pub struct Mirror {
    mode_type: String,
    n_mode: usize,
}
impl Default for Mirror {
    fn default() -> Self {
        Mirror {
            mode_type: String::new(),
            n_mode: 0,
        }
    }
}
pub enum CeoElement {
    Gmt {
        m1: Mirror,
        m2: Mirror,
    },
    Source {
        size: usize,
        pupil_size: f64,
        pupil_sampling: usize,
        band: String,
        zenith: Vec<f32>,
        azimuth: Vec<f32>,
        magnitude: Vec<f32>,
    },
    Shackhartmann {
        n_sensor: usize,
        n_side_lenslet: usize,
        n_px_lenslet: usize,
        d: f64,
        n_px_framelet: usize,
        n_px_imagelet: Option<usize>,
        osf: Option<usize>,
    },
    PSSn {
        r0_at_zenith: f64,
        oscale: f64,
        zenith_angle: f64,
    }
}
pub struct CEO<T> {
    element: CeoElement,
    t: std::marker::PhantomData<T>,
}
pub mod element {
    pub enum ShackHartmann {
        GEOMETRIC,
        DIFFRACTIVE,
    }
    pub struct GMT {}
    pub struct SOURCE {}
    pub struct SHACKHARTMANN {}
    pub struct PSSN {}
}

pub trait Conversion<T> {
    fn from_arcmin(self) -> T;
    fn from_arcsec(self) -> T;
    fn from_mas(self) -> T;
    fn to_arcmin(self) -> T;
    fn to_arcsec(self) -> T;
    fn to_mas(self) -> T;
}
impl Conversion<f64> for f64 {
    /// Converts angle in arcminute to radian
    fn from_arcmin(self) -> f64 {
        self.to_radians() / 60.
    }
    /// Converts angle in arcsecond to radian
    fn from_arcsec(self) -> f64 {
        self.from_arcmin() / 60.
    }
    /// Converts angle in milli-arcsecond to radian
    fn from_mas(self) -> f64 {
        self.from_arcsec() * 1e-3
    }
    /// Converts angle in radian to arcminute
    fn to_arcmin(self) -> f64 {
        60.0 * self.to_degrees()
    }
    /// Converts angle in radian to arcsecond
    fn to_arcsec(self) -> f64 {
        60.0 * self.to_arcmin()
    }
    /// Converts angle in radian to mill-arcsecond
    fn to_mas(self) -> f64 {
        1e3 * self.to_arcsec()
    }
}
impl Conversion<f32> for f32 {
    /// Converts angle in arcminute to radian
    fn from_arcmin(self) -> f32 {
        self.to_radians() / 60.
    }
    /// Converts angle in arcsecond to radian
    fn from_arcsec(self) -> f32 {
        self.from_arcmin() / 60.
    }
    /// Converts angle in milli-arcsecond to radian
    fn from_mas(self) -> f32 {
        self.from_arcsec() * 1e-3
    }
    /// Converts angle in radian to arcminute
    fn to_arcmin(self) -> f32 {
        60.0 * self.to_degrees()
    }
    /// Converts angle in radian to arcsecond
    fn to_arcsec(self) -> f32 {
        60.0 * self.to_arcmin()
    }
    /// Converts angle in radian to mill-arcsecond
    fn to_mas(self) -> f32 {
        1e3 * self.to_arcsec()
    }
}

pub fn set_gpu(id: i32) {
    unsafe {
        set_device(id);
    }
}

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
    pub fn filter(&mut self, f: &mut Cu<f32>) -> &mut Self {
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
