use std::mem;

use super::ceo_bindings::{dev2host,imaging};
use super::Propagation;
use super::Source;

pub struct Imaging {
    _c_: imaging,
}
impl Imaging {
    pub fn new() -> Imaging {
        Imaging {
            _c_: unsafe { mem::zeroed() },
        }
    }
    pub fn build(
        &mut self,
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        dft_osf: i32,
        n_px_imagelet: i32,
        binning: i32,
    ) -> &mut Self {
        unsafe {
            self._c_.setup3(
                n_px_lenslet,
                n_side_lenslet,
                dft_osf,
                n_px_imagelet,
                binning,
                n_sensor,
            );
        }
        self
    }
    pub fn __ceo__(&self) -> &imaging {
        &self._c_
    }
    pub fn frame_transfer(&mut self, frame: &mut Vec<f32>) -> &mut Self {
        unsafe{
            dev2host(frame.as_mut_ptr(),
                     self._c_.d__frame,
                     self.resolution()*self.resolution()*self._c_.N_SOURCE);
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe { self._c_.reset(); }
        self
    }
    pub fn resolution(&self) -> i32 {
        self._c_.N_PX_CAMERA*self._c_.N_SIDE_LENSLET
    }
    pub fn n_frame(&self) -> u32 {
        self._c_.N_FRAME as u32
    }
}
impl Drop for Imaging {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Imaging {
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
