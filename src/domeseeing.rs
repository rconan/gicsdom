use ensure::Present;
use ndarray::Array1;
use ndarray_npy::NpzReader;
use s3_sync::{Bucket, Object, Region, S3};
use std::f64;
use std::io::{Cursor, Read};
use std::time::Instant;

pub struct DomeSeeing {
    bucket: Present<Bucket>,
    zenith: u32,
    azimuth: u32,
    enclosure: String,
    wind_speed: u32,
    s3: S3,
    keys: Vec<String>,
    time: Vec<f64>,
    inc: f64,
    step: f64,
    buffer: Vec<Vec<f64>>,
}

impl DomeSeeing {
    pub fn new(
        zenith: u32,
        azimuth: u32,
        enclosure: &str,
        wind_speed: u32,
        rate: Option<f64>,
    ) -> DomeSeeing {
        DomeSeeing {
            bucket: Present(Bucket::from_name(String::from("gmto.starccm"))),
            zenith,
            azimuth,
            enclosure: enclosure.to_owned(),
            wind_speed,
            inc: 5.0 / rate.or(Some(5.0)).unwrap(),
            step: 0.0,
            s3: S3::new(Region::UsEast2),
            keys: Vec::new(),
            time: Vec::new(),
            buffer: Vec::new(),
        }
    }
    pub fn list(&mut self) -> &mut Self {
        let prefix = format!(
            "Baseline2020/b2019_{}z_{}az_{}_{}ms/OPDData_OPD_Data_",
            self.zenith, self.azimuth, self.enclosure, self.wind_speed
        );
        //println!("Scanning {}*.npz:", prefix);
        let o = self.s3.list_objects(&self.bucket, prefix);
        self.keys = o
            .map(|x| x.unwrap().key().to_owned())
            .filter(|x| x.ends_with(".npz"))
            .collect();
        let n = self.keys.len();
        //println!(" * # of keys: {}", n);
        let f = |s: &str| -> f64 {
            s.trim_end_matches(".npz")
                .split("/")
                .last()
                .unwrap()
                .split("_")
                .last()
                .unwrap()
                .parse::<f64>()
                .unwrap()
        };
        self.time = self.keys.iter().map(|x| f(&x) - f(&self.keys[0])).collect();
        //println!(" * time range: [{};{}]", self.time[0], self.time[n - 1]);
        self
    }
    pub fn load_at(&self, t: f64) -> Vec<f64> {
        print!("Loading ");
        let now = Instant::now();
        let q = self
            .time
            .iter()
            .map(|x| (x - t).abs())
            .fold(f64::INFINITY, |a, x| f64::min(a, x));
        let pos = self.time.iter().position(|&x| (x - t).abs() == q).unwrap();
        let key = self.keys[pos].clone();
        print!(" {}", key);
        let object = Object::from_key(&self.bucket, key);
        let mut body = Vec::new();
        self.s3
            .get_body(&Present(object))
            .expect("object body")
            .read_to_end(&mut body)
            .unwrap();
        let r = Cursor::new(body);
        let mut npz = NpzReader::new(r).unwrap();
        let opd: Array1<f64> = npz.by_name("opd.npy").unwrap();
        println!(" in {}ms", now.elapsed().as_millis());
        opd.to_vec()
    }
    pub fn load(&self, key: &str) -> Vec<f64> {
        //print!("Loading  {}", key);
        let now = Instant::now();
        let object = Object::from_key(&self.bucket, key.to_owned());
        let mut body = Vec::new();
        self.s3
            .get_body(&Present(object))
            .expect("object body")
            .read_to_end(&mut body)
            .unwrap();
        let r = Cursor::new(body);
        let mut npz = NpzReader::new(r).unwrap();
        let opd: Array1<f64> = npz.by_name("opd.npy").unwrap();
        //println!(" in {}ms", now.elapsed().as_millis());
        opd.to_vec()
    }
}
impl Iterator for DomeSeeing {
    type Item = Vec<f64>;
    fn next(&mut self) -> Option<Self::Item> {
        let n = (self.keys.len() - 1) as f64;
        if self.step > n {
            self.inc = -self.inc;
            self.step = n + self.inc;
        }
        if self.step < 0.0 {
            self.inc = -self.inc;
            self.step = self.inc;
        }
        let buffer_pos = (self.step / self.inc.abs()).trunc() as usize;
        if (buffer_pos + 1) > self.buffer.len() {
            let pos = self.step.trunc() as usize;
            self.buffer.push(self.load(&self.keys[pos]));
        }
        /*
        println!(
            "inc: {} ; step: {} ; buf.: [len.: {} ; pos.: {}]",
            self.inc,
            self.step,
            self.buffer.len(),
            buffer_pos
        );
        */
        self.step += self.inc;
        Some(self.buffer[buffer_pos].clone())
    }
}
