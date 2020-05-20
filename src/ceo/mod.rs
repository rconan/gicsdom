use std::{f32, f64, fmt, mem};

pub mod atmosphere;
pub mod calibrations;
pub mod centroiding;
pub mod ceo_bindings;
pub mod gmt;
pub mod imaging;
pub mod shackhartmann;
pub mod source;

pub use self::atmosphere::Atmosphere;
pub use self::calibrations::Calibration;
pub use self::centroiding::Centroiding;
pub use self::gmt::{Gmt, GmtState};
pub use self::imaging::{Imaging, LensletArray};
pub use self::shackhartmann::{GeometricShackHartmann, ShackHartmann};
pub use self::source::Propagation;
pub use self::source::Source;
pub use ceo_bindings::{gpu_double, gpu_float, pssn, set_device};

pub trait Conversion {
    fn from_arcmin(self) -> f64;
    fn from_arcsec(self) -> f64;
    fn from_mas(self) -> f64;
    fn to_arcmin(self) -> f64;
    fn to_arcsec(self) -> f64;
    fn to_mas(self) -> f64;
}
impl Conversion for f64 {
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

pub fn set_gpu(id: i32) {
    unsafe {
        set_device(id);
    }
}

pub struct CuFloat {
    _c_: gpu_float,
}
impl CuFloat {
    pub fn new() -> CuFloat {
        CuFloat {
            _c_: unsafe { mem::zeroed() },
        }
    }
    pub fn malloc(&mut self, len: i32) -> &mut Self {
        unsafe {
            self._c_.setup1(len);
            self._c_.dev_malloc();
        }
        self
    }
    pub fn up(&mut self, host_data: &mut Vec<f32>) -> &mut Self {
        unsafe {
            self._c_.host_data = host_data.as_mut_ptr();
            self._c_.host2dev();
        }
        self
    }
    pub fn as_mut_ptr(&mut self) -> *mut f32 {
        self._c_.dev_data
    }
}
impl Drop for CuFloat {
    fn drop(&mut self) {
        unsafe { self._c_.free_dev() }
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
/// let mut pssn = ceo::PSSn::new(0.15,25.0);
/// pssn.build(&mut src);
/// println!("PSSn: {:?}",pssn.reset(&mut src).estimates);
/// ```
///
pub struct PSSn {
    _c_: pssn,
    r_not: f64,
    l_not: f64,
    /// PSSn estimates
    pub estimates: Vec<f32>,
}
impl PSSn {
    /// Creates a new `PSSm`
    pub fn new(r_not: f64, l_not: f64) -> PSSn {
        PSSn {
            _c_: unsafe { mem::zeroed() },
            r_not,
            l_not,
            estimates: vec![],
        }
    }
    /// Initializes PSSn atmosphere and telescope transfer function from a `Source` object
    pub fn build(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_
                .setup(&mut src._c_, self.r_not as f32, self.l_not as f32);
        }
        self.estimates = vec![0.0; self._c_.N as usize];
        self
    }
    /// Estimates the `PSSn` values and resets the `Source` optical transfer function to its initial value
    pub fn reset(&mut self, src: &mut Source) -> &mut Self {
        self.peek(src);
        self._c_.N_O = 0;
        self
    }
    /// Estimates the `PSSn` values
    pub fn peek(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.otf(&mut src._c_);
            self._c_.eval1(self.estimates.as_mut_ptr())
        }
        self
    }
    /// Integrates the `Source` optical transfer function
    pub fn accumulate(&mut self, src: &mut Source) {
        unsafe {
            self._c_.otf(&mut src._c_);
        }
    }
    /// Computes `PSSn` spatial uniformity
    pub fn spatial_uniformity(&mut self) -> f32 {
        let mut pssn_values = self.estimates.clone();
        pssn_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        100. * ((pssn_values.len() as f32)
            * (*pssn_values.last().unwrap() - *pssn_values.first().unwrap()))
            / pssn_values.iter().sum::<f32>()
    }
}
impl fmt::Display for PSSn {
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
