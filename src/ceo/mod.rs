use std::{f32, mem};

pub mod atmosphere;
pub mod centroiding;
pub mod ceo_bindings;
pub mod gmt;
pub mod imaging;
pub mod shackhartmann;
pub mod source;
pub mod calibrations;

pub use self::atmosphere::Atmosphere;
pub use self::centroiding::Centroiding;
pub use self::gmt::{Gmt, GmtState};
pub use self::imaging::Imaging;
pub use self::shackhartmann::{GeometricShackHartmann, ShackHartmann};
pub use self::source::Propagation;
pub use self::source::Source;
pub use self::calibrations::Calibration;
pub use ceo_bindings::{pssn, set_device, gpu_float, gpu_double};

pub fn set_gpu(id: i32) {
    unsafe {
        set_device(id);
    }
}

pub struct CuFloat {
    _c_: gpu_float
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
        unsafe {
            self._c_.free_dev()
        }
    }
}

pub struct PSSn {
    _c_: pssn,
    r_not: f64,
    l_not: f64,
    pub estimates: Vec<f32>,
}
impl PSSn {
    pub fn new(r_not: f64, l_not: f64) -> PSSn {
        PSSn {
            _c_: unsafe { mem::zeroed() },
            r_not,
            l_not,
            estimates: vec![],
        }
    }
    pub fn build(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_
                .setup(&mut src._c_, self.r_not as f32, self.l_not as f32);
        }
        self.estimates = vec![0.0; self._c_.N as usize];
        self
    }
    pub fn reset(&mut self, src: &mut Source) -> &mut Self {
        self.peek(src);
        self._c_.N_O = 0;
        self
    }
    pub fn peek(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.otf(&mut src._c_);
            self._c_.eval1(self.estimates.as_mut_ptr())
        }
        self
    }
    pub fn accumulate(&mut self, src: &mut Source) {
        unsafe {
            self._c_.otf(&mut src._c_);
        }
    }
    pub fn spatial_uniformity(&mut self) -> f32 {
        let mut pssn_values = self.estimates.clone();
        pssn_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        100. * ((pssn_values.len() as f32)
            * (*pssn_values.last().unwrap() - *pssn_values.first().unwrap()))
            / pssn_values.iter().sum::<f32>()
    }
}
