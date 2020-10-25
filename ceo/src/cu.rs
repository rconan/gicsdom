use super::ceo_bindings::{gpu_double, gpu_float};
use core::ops::Mul;
use std::mem;

pub type Single = gpu_float;
pub type Double = gpu_double;

pub trait CuType {
    fn build(&mut self, size: i32);
    fn malloc(&mut self);
    fn assign_f32(&mut self, _ptr: *mut f32) {}
    fn assign_f64(&mut self, _ptr: *mut f64) {}
    fn drop(&mut self);
}
impl CuType for Double {
    fn build(&mut self, size: i32) {
        unsafe {
            self.setup1(size);
        }
    }
    fn malloc(&mut self) {
        unsafe {
            self.dev_malloc();
        }
    }
    fn assign_f64(&mut self, ptr: *mut f64) {
        self.dev_data = ptr;
    }
    fn drop(&mut self) {
        unsafe {
            self.free_dev();
        }
    }
}
impl CuType for Single {
    fn build(&mut self, size: i32) {
        unsafe {
            self.setup1(size);
        }
    }
    fn malloc(&mut self) {
        unsafe {
            self.dev_malloc();
        }
    }
    fn assign_f32(&mut self, ptr: *mut f32) {
        self.dev_data = ptr;
    }
    fn drop(&mut self) {
        unsafe {
            self.free_dev();
        }
    }
}
pub struct Cu<T: CuType> {
    _c_: T,
    m: usize,
    n: usize,
    dev_alloc: bool,
    marker: std::marker::PhantomData<T>,
}
impl<T: CuType> Cu<T> {
    pub fn new() -> Cu<T> {
        Cu {
            _c_: unsafe { mem::zeroed() },
            m: 0,
            n: 0,
            dev_alloc: false,
            marker: std::marker::PhantomData,
        }
    }
    pub fn array(m: usize, n: usize) -> Cu<T> {
        Cu {
            _c_: unsafe { mem::zeroed() },
            m: m,
            n: n,
            dev_alloc: false,
            marker: std::marker::PhantomData,
        }
    }
    pub fn vector(m: usize) -> Cu<T> {
        Cu {
            _c_: unsafe { mem::zeroed() },
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
    pub fn malloc(&mut self) -> &mut Self {
        let s = self.size();
        self._c_.build(s as i32);
        self._c_.malloc();
        self.dev_alloc = true;
        self
    }
}
impl Cu<Single> {
    pub fn from_ptr(&mut self, ptr: *mut f32) {
        let s = self.size();
        self._c_.build(s as i32);
        self._c_.assign_f32(ptr);
    }
    pub fn as_ptr(&mut self) -> *mut f32 {
        self._c_.dev_data
    }
    pub fn as_mut_ptr(&mut self) -> *mut f32 {
        self._c_.dev_data
    }
    pub fn to_dev(&mut self, host_data: &mut [f32]) -> &mut Self {
        if !self.dev_alloc {
            self.malloc();
        }
        unsafe {
            self._c_.host_data = host_data.as_mut_ptr();
            self._c_.host2dev();
        }
        self
    }
    pub fn from_dev(&mut self) -> Vec<f32> {
        let mut v = vec![0f32; self.size()];
        unsafe {
            self._c_.host_data = v.as_mut_ptr();
            self._c_.dev2host();
        }
        v
    }
    pub fn mv(&mut self, x: &mut Cu<Single>) -> Cu<Single> {
        assert_eq!(x.n_col(), 1, "x must be a vector (n_col=n=1)!");
        assert_eq!(
            x.size(),
            self.n,
            "the number of columns ({}) do not match the length ({})of x",
            self.n,
            x.size()
        );
        let mut y = Cu::<Single>::vector(self.m);
        y.malloc();
        unsafe {
            self._c_.mv(&mut y._c_, &mut x._c_);
        }
        y
    }
    pub fn qr(&mut self) -> &mut Self {
        unsafe {
            self._c_.qr(self.m as i32);
        }
        self
    }
    pub fn qr_solve(&mut self, b: &mut Cu<Single>) -> Cu<Single> {
        let mut x = Cu::<Single>::vector(self.n);
        x.malloc();
        unsafe {
            self._c_.qr_solve(&mut x._c_, &mut b._c_);
        }
        x
    }
    pub fn qr_solve_as_ptr(&mut self, x: &mut Cu<Single>, b: &mut Cu<Single>) {
        unsafe {
            self._c_.qr_solve(&mut x._c_, &mut b._c_);
        }
    }
}
impl Mul for Cu<Single> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut q = rhs;
        let mut s = self;
        s.mv(&mut q)
    }
}
impl Clone for Cu<Single> {
    fn clone(&self) -> Self {
        let mut s = self._c_;
        let mut other = Cu::<Single>::array(self.m, self.n);
        other.malloc();
        unsafe {
            s.dev2dev(&mut other._c_);
        }
        other
    }
}
impl<T: CuType> Drop for Cu<T> {
    fn drop(&mut self) {
        if self.dev_alloc {
            self._c_.drop();
        }
        self.dev_alloc = false;
    }
}
impl From<Vec<f32>> for Cu<Single> {
    fn from(item: Vec<f32>) -> Self {
        let mut this = Cu::vector(item.len());
        this.to_dev(&mut item.clone());
        this
    }
}

impl From<Cu<Single>> for Vec<f32> {
    fn from(item: Cu<Single>) -> Self {
        let mut q = item;
        q.from_dev()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cu_todev() {
        let mut d = Cu::<Single>::vector(5);
        let mut v = vec![1f32; 5];
        d.malloc();
        d.to_dev(&mut v.as_mut_slice());
    }

    #[test]
    fn cu_fromdev() {
        let mut d = Cu::<Single>::vector(5);
        let mut v = vec![1f32; 5];
        d.malloc().to_dev(&mut v.as_mut_slice());
        let w = d.from_dev();
        println!("w: {:?}", w);
    }

    #[test]
    fn cu_toptr() {
        let mut d = Cu::<Single>::vector(5);
        let mut v = vec![4.321f32; 5];
        d.malloc().to_dev(&mut v.as_mut_slice());
        let w = d.from_dev();
        println!("w: {:?}", w);
        let mut dd = Cu::<Single>::vector(5);
        dd.from_ptr(d.as_ptr());
        let w = dd.from_dev();
        println!("w: {:?}", w);
    }

    #[test]
    fn cu_from_into() {
        let d_v = Cu::<Single>::from(vec![1f32; 7]);
        let u: Vec<f32> = d_v.into();
        println!("u: {:?}", u);
        assert_eq!(u, vec![1f32; 7]);
    }

    #[test]
    fn cu_into_from() {
        let v = vec![1f32; 7];
        let d_v: Cu<Single> = v.into();
        let u = Vec::<f32>::from(d_v);
        assert_eq!(u, vec![1f32; 7]);
    }
}
