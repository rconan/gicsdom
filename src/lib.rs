extern crate hifitime;

use crate::ceo::Propagation;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::f32;
use std::f64;
use std::time::Instant;

pub mod ceo;
pub mod astrotools;
pub mod agws;

use astrotools::Observation;

const DEG2RAD: f64 = f64::consts::PI / 180.0;
const RAD2ARCMIN: f64 = 180. * 60. / f64::consts::PI;
const GMT_LAT: f64 = -29.049;
const GMT_LONG: f64 = -70.682;

pub trait AngleConversion {
    fn deg2rad(self) -> f64;
    fn rad2deg(self) -> f64;
}
impl AngleConversion for f64 {
    fn deg2rad(self) -> f64 {
        self * DEG2RAD
    }
    fn rad2deg(self) -> f64 {
        self / DEG2RAD
    }
}

pub struct Rotation {
    pub o: f64,
    pub axis: u8,
}
impl Rotation {
    pub fn new(o: f64, axis: u8) -> Rotation {
        Rotation { o, axis }
    }
    pub fn apply(&self, v: Vec<f64>) -> Vec<f64> {
        let (s, c) = self.o.sin_cos();
        let x = v[0];
        let y = v[1];
        let z = v[2];
        let rot_v: Vec<f64> = match self.axis {
            0 => vec![x, c * y + s * z, -s * y + c * z],
            1 => vec![c * x - s * z, y, s * x + c * z],
            2 => vec![c * x + s * y, -s * x + c * y, z],
            _ => {
                println!("Rotation axis must be one of: 0 for x, 1 for y or 2 for z!");
                vec![0.0; 3]
            }
        };
        rot_v
    }
}


