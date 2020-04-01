use ndarray::{Array1};
use ndarray_npy::NpzReader;
use std::fs::File;
use s3_sync::{S3, Region, Object, Bucket};
use std::io::{Cursor,Read};
use ensure::Present;
use std::time::Instant;

fn main() {
    let mut npz = NpzReader::new(File::open("/home/rconan/Downloads/OPDData_OPD_Data_5.000000e+02.npz").unwrap()).unwrap();
    let names = npz.names().unwrap();
    println!("{:?}",names);
    let opd: Array1<f64> = npz.by_name("opd.npy").unwrap();
    println!("{:?}",opd.shape());
    let iter = opd.iter().filter(|x| !x.is_nan());
    let n = iter.clone().fold(0,|s,x| s + 1 ) as f64;
    let m = iter.clone().fold(0.0,|s,x| s + x )/n;
    let s2 = iter.clone().fold(0.0,|s,x| s + (x-m).powi(2) )/n;
    println!("n={} ; m={} ; s={}",n,m*1e6,1e6*s2.sqrt());

    let s3 = S3::new(Region::UsEast2);
    let bucket = Present(Bucket::from_name("gmto.starccm".to_owned()));
    let object = Object::from_key(&bucket, "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_5.000000e+02.npz".to_owned());
    let mut body = Vec::new();

    let now = Instant::now();
    s3.get_body(&Present(object)).expect("object body").read_to_end(&mut body).unwrap();
    let r = Cursor::new(body);
    let mut npz = NpzReader::new(r).unwrap();
    println!("Loading in {}s", now.elapsed().as_secs());

    let names = npz.names().unwrap();
    println!("{:?}",names);
    let opd: Array1<f64> = npz.by_name("opd.npy").unwrap();
    println!("{:?}",opd.shape());
    let iter = opd.iter().filter(|x| !x.is_nan());
    let n = iter.clone().fold(0,|s,x| s + 1 ) as f64;
    let m = iter.clone().fold(0.0,|s,x| s + x )/n;
    let s2 = iter.clone().fold(0.0,|s,x| s + (x-m).powi(2) )/n;
    println!("n={} ; m={} ; s={}",n,m*1e6,1e6*s2.sqrt());
}
