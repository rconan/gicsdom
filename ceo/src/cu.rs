use super::ceo_bindings::{gpu_double, gpu_float};
use std::mem;
use core::ops::Mul;

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
    pub fn n_row(&mut self) -> usize {
        self.m
    }
    pub fn n_col(&mut self) -> usize {
        self.n
    }
}
impl Cu<f32> {
    pub fn malloc(&mut self) -> &mut Self {
        let s = self.size();
        unsafe {
            self._c_f32.setup1(s as i32);
            self._c_f32.dev_malloc();
        }
        self.dev_alloc = true;
        self
    }
    pub fn from_ptr(&mut self, ptr: *mut f32) {
        let s = self.size();
        unsafe {
            self._c_f32.setup1(s as i32);
            self._c_f32.dev_data = ptr;
        }
    }
    pub fn as_ptr(&mut self) -> *mut f32 {
        self._c_f32.dev_data
    }
    pub fn as_mut_ptr(&mut self) -> *mut f32 {
        self._c_f32.dev_data
    }
    pub fn to_dev(&mut self, host_data: &mut [f32]) -> &mut Self {
        if !self.dev_alloc {
            self.malloc();
        }
        unsafe {
            self._c_f32.host_data = host_data.as_mut_ptr();
            self._c_f32.host2dev();
        }
        self
    }
    pub fn from_dev(&mut self) -> Vec<f32> {
        let mut v = vec![0f32; self.size()];
        unsafe {
            self._c_f32.host_data = v.as_mut_ptr();
            self._c_f32.dev2host();
        }
        v
    }
    pub fn mv(&mut self, x: &mut Cu<f32>) -> Cu<f32> {
        assert_eq!(x.n_col(), 1, "x must be a vector (n_col=n=1)!");
        assert_eq!(
            x.size(),
            self.n,
            "the number of columns ({}) do not match the length ({})of x",
            self.n,
            x.size()
        );
        let mut y = Cu::<f32>::vector(self.m);
        y.malloc();
        unsafe {
            self._c_f32.mv(&mut y._c_f32, &mut x._c_f32);
        }
        y
    }
    pub fn qr(&mut self) -> &mut Self {
        unsafe {
            self._c_f32.qr(self.m as i32);
        }
        self
    }
    pub fn qr_solve(&mut self, b: &mut Cu<f32>) -> Cu<f32> {
        let mut x = Cu::<f32>::vector(self.n);
        x.malloc();
        unsafe {
            self._c_f32.qr_solve(&mut x._c_f32, &mut b._c_f32);
        }
        x
    }
    pub fn qr_solve_as_ptr(&mut self, x: &mut Cu<f32>, b: &mut Cu<f32>) {
        unsafe {
            self._c_f32.qr_solve(&mut x._c_f32, &mut b._c_f32);
        }
    }
}
impl Mul for Cu<f32> {
    // The multiplication of rational numbers is a closed operation.
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut q = rhs;
        let mut s = self;
        s.mv(&mut q)
    }
}
impl Clone for Cu<f32> {
    fn clone(&self) -> Self {
        let mut s  = self._c_f32;
        let mut other = Cu::<f32>::array(self.m, self.n);
        other.malloc();
        unsafe {
            s.dev2dev(&mut other._c_f32);
        }
        other
    }
}
impl<T> Drop for Cu<T> {
    fn drop(&mut self) {
        if self.dev_alloc {
            unsafe { self._c_f32.free_dev() }
        }
        self.dev_alloc = false;
    }
}
impl From<Vec<f32>> for Cu<f32> {
    fn from(item: Vec<f32>) -> Self {
        let mut this = Cu::vector(item.len());
        this.to_dev(&mut item.clone());
        this
    }
}

impl From<Cu<f32>> for Vec<f32> {
    fn from(item: Cu<f32>) -> Self {
        let mut q = item;
        q.from_dev()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cu_todev() {
        let mut d = Cu::<f32>::vector(5);
        let mut v = vec![1f32; 5];
        d.malloc();
        d.to_dev(&mut v.as_mut_slice());
    }

    #[test]
    fn cu_fromdev() {
        let mut d = Cu::<f32>::vector(5);
        let mut v = vec![1f32; 5];
        d.malloc().to_dev(&mut v.as_mut_slice());
        let w = d.from_dev();
        println!("w: {:?}", w);
    }

    #[test]
    fn cu_toptr() {
        let mut d = Cu::<f32>::vector(5);
        let mut v = vec![4.321f32; 5];
        d.malloc().to_dev(&mut v.as_mut_slice());
        let w = d.from_dev();
        println!("w: {:?}", w);
        let mut dd = Cu::<f32>::vector(5);
        dd.from_ptr(d.as_ptr());
        let w = dd.from_dev();
        println!("w: {:?}", w);
    }

    #[test]
    fn cu_from_into() {
        let d_v = Cu::<f32>::from(vec![1f32; 7]);
        let u: Vec<f32> = d_v.into();
        println!("u: {:?}", u);
        assert_eq!(u, vec![1f32; 7]);
    }

    #[test]
    fn cu_into_from() {
        let v = vec![1f32; 7];
        let d_v: Cu<f32> = v.into();
        let u = Vec::<f32>::from(d_v);
        assert_eq!(u, vec![1f32; 7]);
    }
}
