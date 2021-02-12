use super::ceo_bindings::{gpu_double, gpu_float};
use core::ops::{AddAssign, Mul, SubAssign};
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
    /// number of rows
    pub n_rows: usize,
    /// number of columns
    pub n_cols: usize,
    dev_alloc: bool,
}
impl<T: CuType> Cu<T> {
    /// Creates an empty CUDA array
    pub fn new() -> Cu<T> {
        Cu {
            _c_: unsafe { mem::zeroed() },
            n_rows: 0,
            n_cols: 0,
            dev_alloc: false,
        }
    }
    pub fn array(n_rows: usize, n_cols: usize) -> Cu<T> {
        Cu {
            _c_: unsafe { mem::zeroed() },
            n_rows: n_rows,
            n_cols: n_cols,
            dev_alloc: false,
        }
    }
    pub fn vector(n_rows: usize) -> Cu<T> {
        Cu {
            _c_: unsafe { mem::zeroed() },
            n_rows: n_rows,
            n_cols: 1,
            dev_alloc: false,
        }
    }
    pub fn size(&self) -> usize {
        self.n_rows * self.n_cols
    }
    pub fn n_rows(&self) -> usize {
        self.n_rows
    }
    pub fn n_cols(&self) -> usize {
        self.n_cols
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
    pub fn mv(&self, x: &Cu<Single>) -> Cu<Single> {
        assert_eq!(x.n_cols(), 1, "x must be a vector (n_col=n=1)!");
        assert_eq!(
            x.size(),
            self.n_cols,
            "the number of columns ({}) do not match the length ({})of x",
            self.n_cols,
            x.size()
        );
        let mut y = Cu::<Single>::vector(self.n_rows);
        y.malloc();
        unsafe {
            let mut s = self._c_;
            s.mv(&mut y._c_, &x._c_);
        }
        y
    }
    pub fn mv_unchecked(&self, y: &mut Cu<Single>, x: &Cu<Single>) {
        unsafe {
            let mut s = self._c_;
            s.mv(&mut y._c_, &x._c_);
        }
    }
    pub fn qr(&mut self) -> &mut Self {
        unsafe {
            self._c_.qr(self.n_rows as i32);
        }
        self
    }
    pub fn qr_solve(&mut self, b: &mut Cu<Single>) -> Cu<Single> {
        let mut x = Cu::<Single>::vector(self.n_cols);
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
impl Mul for &Cu<Single> {
    type Output = Cu<Single>;

    fn mul(self, rhs: Self) -> Cu<Single> {
        self.mv(rhs)
    }
}
impl Mul<f32> for &Cu<Single> {
    type Output = Cu<Single>;

    fn mul(self, rhs: f32) -> Cu<Single> {
        let mut other = self.clone();
        unsafe {
            other._c_.scale(rhs);
        }
        other
    }
}
impl AddAssign for Cu<Single> {
    fn add_assign(&mut self, other: Self) {
        unsafe {
            self._c_.axpy(&other._c_, 1.0);
        }
    }
}
impl SubAssign for Cu<Single> {
    fn sub_assign(&mut self, other: Self) {
        unsafe {
            self._c_.axpy(&other._c_, -1.0);
        }
    }
}
impl Clone for Cu<Single> {
    fn clone(&self) -> Self {
        let mut s = self._c_;
        let mut other = Cu::<Single>::array(self.n_rows, self.n_cols);
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
        let mut this = Cu::<Single>::vector(item.len());
        this.to_dev(&mut item.clone());
        this
    }
}
impl From<Vec<Vec<f32>>> for Cu<Single> {
    fn from(item: Vec<Vec<f32>>) -> Self {
        let n_cols = item.len();
        let n_rows = item[0].len();
        let mut flat_item: Vec<f32> = item.iter().cloned().flatten().collect();
        let mut this = Cu::<Single>::array(n_rows, n_cols);
        this.to_dev(&mut flat_item);
        this
    }
}

impl From<Cu<Single>> for Vec<f32> {
    fn from(item: Cu<Single>) -> Self {
        let mut q = item;
        q.from_dev()
    }
}

impl From<&mut Cu<Single>> for Vec<f32> {
    fn from(item: &mut Cu<Single>) -> Self {
        item.from_dev()
    }
}

impl Cu<Double> {
    pub fn from_dev(&mut self) -> Vec<f64> {
        let mut v = vec![0f64; self.size()];
        unsafe {
            self._c_.host_data = v.as_mut_ptr();
            self._c_.dev2host();
        }
        v
    }
    pub fn to_dev(&mut self, host_data: &mut [f64]) -> &mut Self {
        if !self.dev_alloc {
            self.malloc();
        }
        unsafe {
            self._c_.host_data = host_data.as_mut_ptr();
            self._c_.host2dev();
        }
        self
    }
}
impl From<Vec<f64>> for Cu<Double> {
    fn from(item: Vec<f64>) -> Self {
        let mut this = Cu::<Double>::vector(item.len());
        this.to_dev(&mut item.clone());
        this
    }
}
impl From<Vec<Vec<f64>>> for Cu<Double> {
    fn from(item: Vec<Vec<f64>>) -> Self {
        let n_cols = item.len();
        let n_rows = item[0].len();
        let mut flat_item: Vec<f64> = item.iter().cloned().flatten().collect();
        let mut this = Cu::<Double>::array(n_rows, n_cols);
        this.to_dev(&mut flat_item);
        this
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
