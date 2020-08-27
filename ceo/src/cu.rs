use super::ceo_bindings::{gpu_double, gpu_float};
use std::mem;

pub struct Cu<T> {
    _c_f32: gpu_float,
    _c_f64: gpu_double,
    m: usize,
    n: usize,
    dev_alloc: bool,
    marker: std::marker::PhantomData<T>,
}
impl<T> Cu<T> {
    pub fn new() -> Cu<T> {
        Cu {
            _c_f32: unsafe { mem::zeroed() },
            _c_f64: unsafe { mem::zeroed() },
            m: 0,
            n: 0,
            dev_alloc: false,
            marker: std::marker::PhantomData,
        }
    }
    pub fn array(m: usize, n: usize) -> Cu<T> {
        Cu {
            _c_f32: unsafe { mem::zeroed() },
            _c_f64: unsafe { mem::zeroed() },
            m: m,
            n: n,
            dev_alloc: false,
            marker: std::marker::PhantomData,
        }
    }
    pub fn vector(m: usize) -> Cu<T> {
        Cu {
            _c_f32: unsafe { mem::zeroed() },
            _c_f64: unsafe { mem::zeroed() },
            m: m,
            n: 1,
            dev_alloc: false,
            marker: std::marker::PhantomData,
        }
    }
    pub fn size(&mut self) -> usize {
        self.m * self.n
    }
}
impl Cu<f32> {
    pub fn build(&mut self) -> &mut Self {
        let s = self.size();
        unsafe {
            self._c_f32.setup1(s as i32);
        }
        self
    }
    pub fn malloc(&mut self) -> &mut Self {
        unsafe {
            self._c_f32.dev_malloc();
        }
        self.dev_alloc = true;
        self
    }
    pub fn to_ptr(&mut self, ptr: *mut f32) {
        self._c_f32.dev_data = ptr;
    }
    pub fn as_ptr(&mut self) -> *mut f32 {
        self._c_f32.dev_data
    }
    pub fn to_dev(&mut self, host_data: &mut [f32]) -> Result<&mut Self,String> {
        if self.dev_alloc {
            unsafe {
                self._c_f32.host_data = host_data.as_mut_ptr();
                self._c_f32.host2dev();
            }
            return Ok(self);
        } else {
            return Err("Memory on device has not been allocated!".to_string());
        }
    }
    pub fn from_dev(&mut self) -> Vec<f32> {
        let mut v = vec![0f32;self.size()];
        unsafe {
            self._c_f32.host_data = v.as_mut_ptr();
            self._c_f32.dev2host();
        }
        v
    }
}
impl<T> Drop for Cu<T> {
    fn drop(&mut self) {
        if self.dev_alloc {
            unsafe { self._c_f32.free_dev() }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cu_todev() {
        let mut d = Cu::<f32>::vector(5);
        let mut v = vec![1f32;5];
        d.build().malloc();
        d.to_dev(&mut v.as_mut_slice()).expect("Failed");
    }

    #[test]
    fn cu_fromdev() {
        let mut d = Cu::<f32>::vector(5);
        let mut v = vec![1f32;5];
        d.build().malloc().to_dev(&mut v.as_mut_slice()).expect("Failed");
        let w = d.from_dev();
        println!("w: {:?}",w);
    }

    #[test]
    fn cu_toptr() {
        let mut d = Cu::<f32>::vector(5);
        let mut v = vec![4.321f32;5];
        d.build().malloc().to_dev(&mut v.as_mut_slice()).expect("Failed");
        let w = d.from_dev();
        println!("w: {:?}",w);
        let mut dd = Cu::<f32>::vector(5);
        dd.build().to_ptr(d.as_ptr());
        let w = dd.from_dev();
        println!("w: {:?}",w);
    }

}
