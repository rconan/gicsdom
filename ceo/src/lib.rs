use std::{f32, f64, fmt, mem};

pub mod atmosphere;
//pub mod calibrations;
pub mod centroiding;
pub mod ceo_bindings;
pub mod gmt;
pub mod imaging;
pub mod shackhartmann;
pub mod source;
pub mod cu;

pub use self::atmosphere::Atmosphere;
//pub use self::calibrations::Calibration;
pub use self::centroiding::Centroiding;
pub use self::gmt::Gmt;
pub use self::imaging::{Imaging, LensletArray};
pub use self::shackhartmann::{GeometricShackHartmann, ShackHartmann};
pub use self::source::Propagation;
pub use self::source::Source;
pub use self::cu::Cu;
pub use ceo_bindings::{gpu_double, gpu_float, pssn, set_device, geqrf, ormqr};

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

pub struct CuFloat {
    _c_: gpu_float,
    host_data: Vec<f32>,
}
impl CuFloat {
    pub fn new() -> CuFloat {
        CuFloat {
            _c_: unsafe { mem::zeroed() },
            host_data: vec![]
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
    unsafe {
        geqrf(tau._c_.dev_data, a._c_.dev_data, m, n)
    }
}
pub fn qtb(b: &mut CuFloat, m: i32, a: &mut CuFloat, tau: &mut CuFloat, n: i32) {
    unsafe {
        ormqr(b._c_.dev_data, m, a._c_.dev_data, tau._c_.dev_data, n)
    }
}

/// Wrapper for CEO PSSn
///
/// # Examples
///
/// ```
/// use gicsdom::ceo;
/// let mut src = ceo::Source::new(1,25.5,401);
/// src.build("V",vec![0.0],vec![0.0],vec![0.0]);
/// let mut gmt = ceo::Gmt::new();
/// gmt.build(27,None);
/// src.through(&mut gmt).xpupil();
/// println!("WFE RMS: {:.3}nm",src.wfe_rms_10e(-9)[0]);
/// let mut pssn = ceo::PSSn::new();
/// pssn.build(&mut src);
/// println!("PSSn: {:?}",pssn.reset(&mut src).estimates);
/// ```
///
// NEW PSSN
pub struct TelescopeError;
pub struct AtmosphereTelescopeError;
pub struct GPSSn<S> {
    _c_: pssn,
    pub r0_at_zenith: f32,
    pub oscale: f32,
    pub zenith_angle: f32,
    /// GPSSn estimates
    pub estimates: Vec<f32>,
    mode: std::marker::PhantomData<S>,
}
impl<S> GPSSn<S> {
    /// Creates a new `GPSSn` with r0=16cm at zenith, L0=25m a zenith distance of 30 degrees
    pub fn new() -> GPSSn<S> {
        GPSSn {
            _c_: unsafe { mem::zeroed() },
            r0_at_zenith: 0.16,
            oscale: 25.0,
            zenith_angle: 30_f32.to_radians(),
            estimates: vec![],
            mode: std::marker::PhantomData,
        }
    }
    /// Initializes GPSSn atmosphere and telescope transfer function from a `Source` object
    pub fn build(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.setup(&mut src._c_, self.r0(), self.oscale);
        }
        self.estimates = vec![0.0; self._c_.N as usize];
        self
    }
    /// Integrates the `Source` optical transfer function
    pub fn accumulate(&mut self, src: &mut Source) {
        unsafe {
            self._c_.otf(&mut src._c_);
        }
    }
    pub fn reset(&mut self) -> &mut Self {
        self._c_.N_O = 0;
        self
    }
    /// Computes `GPSSn` spatial uniformity
    pub fn spatial_uniformity(&mut self) -> f32 {
        let mut pssn_values = self.estimates.clone();
        pssn_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        100. * ((pssn_values.len() as f32)
            * (*pssn_values.last().unwrap() - *pssn_values.first().unwrap()))
            / pssn_values.iter().sum::<f32>()
    }
    pub fn r0(&self) -> f32 {
        (self.r0_at_zenith.powf(-5_f32 / 3_f32) / self.zenith_angle.cos()).powf(-3_f32 / 5_f32)
    }
}
impl GPSSn<TelescopeError> {
    /// Estimates the `GPSSn` values
    pub fn peek(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.eval1(self.estimates.as_mut_ptr())
        }
        self
    }
}
impl GPSSn<AtmosphereTelescopeError> {
    /// Estimates the `GPSSn` values
    pub fn peek(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.oeval1(self.estimates.as_mut_ptr())
        }
        self
    }
}
impl<S> Drop for GPSSn<S> {
    /// Frees CEO memory before dropping `GPSSn`
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl<S> fmt::Display for GPSSn<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}]",
            self.estimates
                .iter()
                .map(|x| format!("{:.4}", x))
                .collect::<Vec<String>>()
                .as_slice()
                .join(",")
        )
    }
}

pub type PSSn = GPSSn<TelescopeError>;

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
                1f32, 4f32, 2f32, 1f32, 2f32, 5f32, 1f32, 1f32, 3f32, 6f32, 1f32, 1f32
            ])
            .qr(4);
        println!("a: {:?}",a.from());
        let mut b = CuFloat::new();
        b.malloc(4).into(&mut vec![6f32,15f32,4f32,3f32]);
        let mut x = CuFloat::new();
        x.malloc(3);
        a.qr_solve(&mut x, &mut b);
        println!("x: {:?}",x.from());
    }

    #[test]
    fn cufloat_bigqr() {
        let m = 5000;
        let n = 700;
        let p = m*n;
        let mut rng = rand::thread_rng();
        let mut a = CuFloat::new();
        a.malloc(p)
            .into(&mut (0..p).map(|_| rng.gen_range(-100,100) as f32).collect::<Vec<f32>>());
        //println!("a: {:?}",a.from());
        let mut b = CuFloat::new();
        b.malloc(m);
        {
            let mut x = CuFloat::new();
            x.malloc(n).into(&mut vec![1f32;n]);
            a.mv(&mut b, &mut x);
        }
        //println!("b: {:?}",b.from());
        a.qr(m as i32);
        //println!("a: {:?}",a.from());
        let mut x = CuFloat::new();
        x.malloc(n);
        a.qr_solve(&mut x, &mut b);
        let sx = x.from().iter().sum::<f32>();
        println!("sum of x: {:?}",sx);
        assert!((sx-(n as f32)).abs()<1e-6);
    }
}
