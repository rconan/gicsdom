extern crate hifitime;

use hifitime::Epoch;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::f32;
use std::f64;
use std::time::Instant;

pub mod bindings;
pub mod ceo;

const DEG2RAD: f64 = f64::consts::PI / 180.0;

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
    pub mat: [[f64; 3]; 3],
}
impl Rotation {
    pub fn new(o: f64, axis: u8) -> Rotation {
        Rotation {
            o,
            mat: {
                let osc = o.sin_cos();
                let z = 0f64;
                match axis {
                    0 => [[1f64, z, z], [z, osc.1, osc.0], [z, -osc.0, osc.1]],
                    1 => [[osc.1, z, -osc.0], [z, 1f64, z], [osc.0, z, osc.1]],
                    2 => [[osc.1, osc.0, z], [-osc.0, osc.1, z], [z, z, 1f64]],
                    _ => {
                        println!("Rotation axis must be one of: 0 for x, 1 for y or 2 for z!");
                        [[z; 3]; 3]
                    }
                }
            },
        }
    }
    pub fn apply(&self, v: Vec<f64>) -> Vec<f64> {
        self.mat
            .iter()
            .map(|x| {
                x.iter()
                    .zip(v.iter())
                    .map(|y| y.0 * y.1)
                    .fold(0.0, |acc, z| acc + z)
            })
            .collect()
    }
}

pub struct Observation {
    pub datetime: String,
    /// site latitude [radian]
    pub latitude: f64,
    /// site longitude [radian]
    pub longitude: f64,
    //    height: f64,
    epoch: Epoch,
}
impl Observation {
    pub fn new(_datetime: &str, latitude_deg: f64, longitude_deg: f64) -> Observation {
        Observation {
            datetime: _datetime.to_string(),
            latitude: latitude_deg.deg2rad(),
            longitude: longitude_deg.deg2rad(),
            epoch: Epoch::from_gregorian_utc_str(_datetime).unwrap(),
        }
    }
    pub fn latitude_deg(&self) -> f64 {
        self.latitude.rad2deg()
    }
    pub fn longitude_deg(&self) -> f64 {
        self.longitude.rad2deg()
    }
    pub fn j2000_day(&self) -> f64 {
        let jde: f64 = self.epoch.as_jde_tt_days();
        jde - 2451545.0
    }
    pub fn local_sidereal_time(&self) -> f64 {
        let mut lst = 100.46
            + 0.985647 * self.j2000_day()
            + self.longitude_deg()
            + 15.0 * self.decimal_hour();
        if lst < 0_f64 {
            lst += 360.0;
        }
        lst.deg2rad()
    }
    pub fn local_sidereal_time_deg(&self) -> f64 {
        self.local_sidereal_time().rad2deg()
    }
    pub fn decimal_hour(&self) -> f64 {
        let (_y, _m, _d, h, min, sec) = self.epoch.as_gregorian_utc();
        (h as f64) + ((min as f64) + (sec as f64) / 60.0) / 60.0
    }
    pub fn add_seconds(&mut self, value: f64) -> &mut Self {
        let tai_sec = self.epoch.as_tai_seconds() + value;
        self.epoch = Epoch::from_tai_seconds(tai_sec);
        self.datetime = self.epoch.as_gregorian_utc_str();
        self
    }
}
pub struct SkyCoordinates {
    /// right ascension and declination [degree]
    pub radec: (f64, f64),
}
pub fn pol2cart(v: Vec<f64>) -> Vec<f64> {
    vec![v[0].hypot(v[1]),v[1].atan2(v[0])]
}
impl SkyCoordinates {
    pub fn new(ra_deg: f64, dec_deg: f64) -> SkyCoordinates {
        SkyCoordinates {
            radec: (ra_deg.deg2rad(), dec_deg.deg2rad()),
        }
    }
    pub fn radec_deg(&self) -> (f64, f64) {
        (self.radec.0.rad2deg(), self.radec.1.rad2deg())
    }
    pub fn hour_angle(&self, obs: &Observation) -> f64 {
        obs.local_sidereal_time() - self.radec.0
    }
    pub fn alt_az_parallactic(&self, obs: &Observation) -> (f64, f64, f64) {
        let sin_cos_ha = self.hour_angle(obs).sin_cos();
        let sin_alt = self.radec.1.sin() * obs.latitude.sin()
            + self.radec.1.cos() * obs.latitude.cos() * sin_cos_ha.1;
        let alt = sin_alt.asin();
        let sin_az = sin_cos_ha.0;
        let cos_az = sin_cos_ha.1 * obs.latitude.sin() - self.radec.1.tan() * obs.latitude.cos();
        let az = sin_az.atan2(cos_az);
        let cos_p = obs.latitude.tan() * self.radec.1.cos() - self.radec.1.sin() * sin_cos_ha.1;
        let p = sin_az.atan2(cos_p);
        (alt, az, p)
    }
    pub fn altaz(&self, obs: &Observation) -> (f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0, aap.1)
    }
    pub fn altaz_deg(&self, obs: &Observation) -> (f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0.rad2deg(), aap.1.rad2deg())
    }
    pub fn local(&self, obs: &Observation) -> Vec<f64> {
        let altaz = self.altaz(obs);
        let sc_alt = altaz.0.sin_cos();
        let sc_az = altaz.1.sin_cos();
        let xyz = vec![sc_az.1 * sc_alt.1, sc_az.0 * sc_alt.1, sc_alt.0];
        Rotation::new(f64::consts::FRAC_PI_2 - altaz.0, 1)
            .apply(Rotation::new(altaz.1, 2).apply(xyz))
    }
}

pub struct OpticalPathToSH48 {
    pub gmt: ceo::Gmt,
    pub gs: ceo::Source,
    pub sensor: ceo::GeometricShackHartmann,
}
impl OpticalPathToSH48 {
    pub fn new(n_sensor: i32) -> OpticalPathToSH48 {
        let n_side_lenslet = 48;
        let n_px_lenslet = 16;
        //let n_px = n_side_lenslet*16 + 1;
        let pupil_size = 25.5;
        let m1_n_mode = 27;
        let d = pupil_size / n_side_lenslet as f64;
        OpticalPathToSH48 {
            gmt: ceo::Gmt::new(m1_n_mode, None),
            gs: ceo::Source::empty(),
            sensor: ceo::Geometric_ShackHartmann::new(n_sensor, n_side_lenslet, n_px_lenslet, d),
        }
    }
    pub fn build(&mut self, zen: Vec<f32>, azi: Vec<f32>) -> &mut Self {
        self.gmt.build();
        self.sensor.build();
        self.gs = self.sensor.new_guide_stars();
        self.gs.build(
            "V",
            zen,
            azi,
            vec![0.0, 0.0, 0.0],
        );
        self.gs.through(&mut self.gmt);
        self.sensor.calibrate(&mut self.gs, 0.0).unwrap();
        self
    }
    pub fn propagate_src(&mut self) {
        self.gs.through(&mut self.gmt).through(&mut self.sensor);
    }
    pub fn calibrate(&mut self, progress: Option<ProgressBar>) -> Array2<f32> {
        let now = Instant::now();
        let n_rbm: usize = 84;
        let mut c_c_p: Vec<f32> = Vec::new();
        let mut c_c_m: Vec<f32> = Vec::new();
        let n_c: usize =
            (self.sensor.n_side_lenslet.pow(2) as usize) * 2 * self.sensor.n_sensor as usize;
        self.gmt.reset();
        let pb = match progress {
            Some(_pb) => _pb,
            None => ProgressBar::new(7),
        };
        pb.set_style(ProgressStyle::default_bar().template("{bar:40.cyan/blue} {msg}"));
        if n_rbm > 0 {
            for mid in 1..3 {
                //pb.set_message(&format!("M{} RBM", mid));
                //pb.reset();
                for sid in 1..8 {
                    //pb.inc(1);
                    for tr in 1..3 {
                        for a in 0..3 {
                            self.gmt.reset();
                            let mut t_xyz = vec![0.0; 3];
                            let mut r_xyz = vec![0.0; 3];
                            if tr == 1 {
                                t_xyz[a] = 1e-6;
                            }
                            if tr == 2 {
                                r_xyz[a] = 1e-6;
                            }
                            if mid == 1 {
                                self.gmt.set_m1_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            if mid == 2 {
                                self.gmt.set_m2_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            self.propagate_src();
                            self.sensor.process();
                            c_c_p.append(&mut self.sensor.centroids);
                        }
                    }
                }
            }
        }

        //let pb_bm = ProgressBar::new(self.gmt.m1_n_mode as u64);
        //pb_bm.set_style(ProgressStyle::default_bar().template("{bar:40.cyan/blue} {msg}"));
        if self.gmt.m1_n_mode > 0 {
            //pb.set_length(self.gmt.m1_n_mode as u64);
            //println!("Bending modes ...");
            let mut a: Vec<f64> = vec![0.0; 7 * self.gmt.m1_n_mode as usize];
            self.gmt.reset();
            for sid in 0..7 {
                //pb.set_message(&format!("S{} BM", sid + 1));
                //pb.reset();
                for k_a in 0..self.gmt.m1_n_mode {
                    //pb.inc(1);
                    let k = k_a + sid * self.gmt.m1_n_mode;
                    a[k as usize] = 1e-6;
                    self.gmt.set_m1_modes(&mut a);
                    self.propagate_src();
                    self.sensor.process();
                    c_c_p.append(&mut self.sensor.centroids);
                    a[k as usize] = 0.0;
                    self.gmt.set_m1_modes(&mut a);
                }
            }
        }

        let _d_p = Array2::from_shape_vec((c_c_p.len() / n_c, n_c), c_c_p)
            .unwrap()
            .t()
            .to_owned()
            * 1e6;

        //        println!("Rigid body motion ...");
        //pb.set_length(7);
        if n_rbm > 0 {
            for mid in 1..3 {
                //pb.set_message(&format!("M{} RBM", mid));
                //pb.reset();
                for sid in 1..8 {
                    //pb.inc(1);
                    for tr in 1..3 {
                        for a in 0..3 {
                            self.gmt.reset();
                            let mut t_xyz = vec![0.0; 3];
                            let mut r_xyz = vec![0.0; 3];
                            if tr == 1 {
                                t_xyz[a] = -1e-6;
                            }
                            if tr == 2 {
                                r_xyz[a] = -1e-6;
                            }
                            if mid == 1 {
                                self.gmt.set_m1_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            if mid == 2 {
                                self.gmt.set_m2_segment_state(sid, &t_xyz, &r_xyz);
                            }
                            self.propagate_src();
                            self.sensor.process();
                            //let sum2: f32 = c.iter().sum();
                            //println!("Centroids sum: {:.16}", sum2);
                            //println!("Diff. centroids sum: {:.14}", sum2 - sum1);
                            c_c_m.append(&mut self.sensor.centroids);
                        }
                    }
                }
            }
        }

        if self.gmt.m1_n_mode > 0 {
            //pb.set_length(self.gmt.m1_n_mode as u64);
            //println!("Bending modes ...");
            let mut a: Vec<f64> = vec![0.0; 7 * self.gmt.m1_n_mode as usize];
            self.gmt.reset();
            for sid in 0..7 {
                //pb.set_message(&format!("S{} BM", sid + 1));
                //pb.reset();
                for k_a in 0..self.gmt.m1_n_mode {
                    //pb.inc(1);
                    let k = k_a + sid * self.gmt.m1_n_mode;
                    a[k as usize] = -1e-6;
                    self.gmt.set_m1_modes(&mut a);
                    self.propagate_src();
                    self.sensor.process();
                    c_c_m.append(&mut self.sensor.centroids);
                    a[k as usize] = 0.0;
                    self.gmt.set_m1_modes(&mut a);
                }
            }
        }
        //pb.finish();

        //println!("WFS centroids #     : {}", self.sensor.n_centroids);
        //println!("WFS centroids length: {}", self.sensor.centroids.len());
        //println!("CCM length: {}", c_c_m.len());

        let _d_m = Array2::from_shape_vec((c_c_m.len() / n_c, n_c), c_c_m)
            .unwrap()
            .t()
            .to_owned()
            * 1e6;
        (_d_p - _d_m) * 0.5
        //pb.println(&format!(" in {}s", now.elapsed().as_secs()));

        //println!("{:?}", _d);
        //println!("shape={:?}, strides={:?}", _d.shape(), _d.strides());
        //        println!("d sum: {}",_d.into.sum());

    }
}
impl Drop for OpticalPathToSH48 {
    fn drop(&mut self) {
        drop(&mut self.gmt);
        drop(&mut self.gs);
        drop(&mut self.sensor);
    }
}
