use std::{f32, mem};

pub mod atmosphere;
pub mod gmt;
pub mod source;
pub mod shackhartmann;
pub mod ceo_bindings;

pub use self::atmosphere::Atmosphere;
pub use self::gmt::Gmt;
pub use self::gmt::GmtState;
pub use self::source::Propagation;
pub use self::source::Source;
pub use self::shackhartmann::{Geometric_ShackHartmann,GeometricShackHartmann,ShackHartmann};
pub use ceo_bindings::{set_device,pssn};

pub fn set_gpu(id: i32) {
    unsafe {
        set_device(id);
    }
}

pub struct PSSn {
    _c_: pssn,
    r_not: f64,
    l_not: f64,
    zenith: f64,
    pub estimates: Vec<f32>,
}
impl PSSn {
    pub fn new(r_not: f64, l_not: f64, zenith: f64) -> PSSn {
        PSSn {
            _c_: unsafe { mem::zeroed() },
            r_not,
            l_not,
            zenith,
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
