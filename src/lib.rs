extern crate hifitime;

use crate::ceo::Propagation;
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
        //let tai_sec = self.epoch.as_tai_seconds() + value;
        self.epoch.mut_add_secs(value); // = Epoch::from_tai_seconds(tai_sec);
        self.datetime = self.epoch.as_gregorian_utc_str();
        self
    }
}

pub struct SkyCoordinates {
    /// right ascension and declination [degree]
    pub radec: (f64, f64),
}
pub fn pol2cart(v: Vec<f64>) -> Vec<f64> {
    let v_pol: Vec<f64> = vec![v[0].hypot(v[1]), v[1].atan2(v[0])];
    v_pol
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
    pub fn alt_az_parallactic_deg(&self, obs: &Observation) -> (f64, f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0.rad2deg(), aap.1.rad2deg(), aap.2.rad2deg())
    }
    pub fn local(&self, reference: &SkyCoordinates, obs: &Observation) -> Vec<f64> {
        let altaz = self.altaz(obs);
        let sc_alt = altaz.0.sin_cos();
        let sc_az = altaz.1.sin_cos();
        let (x, y, z) = (sc_az.1 * sc_alt.1, sc_az.0 * sc_alt.1, sc_alt.0);
        let ref_altaz = reference.altaz(obs);
        /*Rotation::new(f64::consts::FRAC_PI_2 - ref_altaz.0, 1)
           .apply(Rotation::new(ref_altaz.1, 2).apply(xyz))
        */
        let (s, c) = ref_altaz.1.sin_cos();
        let w: Vec<f64> = vec![c * x + s * y, -s * x + c * y, z];
        let (s, c) = (f64::consts::FRAC_PI_2 - ref_altaz.0).sin_cos();
        let (x, y, z) = (w[0], w[1], w[2]);
        let w: Vec<f64> = vec![c * x - s * z, y, s * x + c * z];
        w
    }
    pub fn local_polar(&self, reference: &SkyCoordinates, obs: &Observation) -> (f32, f32) {
        let u = self.local(reference, obs);
        let v = pol2cart(u);
        (v[0] as f32, v[1] as f32)
    }
    pub fn local_polar_arcmin_deg(
        &self,
        reference: &SkyCoordinates,
        obs: &Observation,
    ) -> (f64, f64) {
        let (r, o) = self.local_polar(reference, obs);
        (r as f64 * RAD2ARCMIN, o as f64 / DEG2RAD)
    }
}

pub struct OpticalPathToSH48 {
    pub obs: Observation,
    pub telescope: SkyCoordinates,
    pub probe: SkyCoordinates,
    pub gmt: ceo::Gmt,
    pub gs: ceo::Source,
    pub sensor: ceo::GeometricShackHartmann,
    pub probe_id: i32,
    alt_az_pa: (f64, f64, f64),
}
impl OpticalPathToSH48 {
    pub fn new(
        datetime: &str,
        telescope_radec_deg: (f64, f64),
        probe_radec_deg: (f64, f64),
        probe_id: i32,
    ) -> OpticalPathToSH48 {
        let n_side_lenslet = 48;
        let n_px_lenslet = 16;
        //let n_px = n_side_lenslet*16 + 1;
        let pupil_size = 25.5;
        let m1_n_mode = 27;
        let d = pupil_size / n_side_lenslet as f64;
        OpticalPathToSH48 {
            obs: Observation::new(datetime, GMT_LAT, GMT_LONG),
            telescope: SkyCoordinates::new(telescope_radec_deg.0, telescope_radec_deg.1),
            probe: SkyCoordinates::new(probe_radec_deg.0, probe_radec_deg.1),
            gmt: ceo::Gmt::new(m1_n_mode, None),
            gs: ceo::Source::empty(),
            sensor: ceo::Geometric_ShackHartmann::new(1, n_side_lenslet, n_px_lenslet, d),
            probe_id,
            alt_az_pa: (0., 0., 0.),
        }
    }
    pub fn build(&mut self, mag: f64) -> &mut Self {
        self.gmt.build();
        self.sensor.build();
        self.gs = self.sensor.new_guide_stars();
        let q = self.probe.local(&self.telescope, &self.obs);
        let z = q[0].hypot(q[1]) as f32;
        let a = q[1].atan2(q[0]) as f32;
        let zen: Vec<f32> = vec![z]; //6.0*f32::consts::PI/180./60.];
        let azi: Vec<f32> = vec![a]; //2. * (self.probe_id as f32) * f32::consts::PI / 3.];
        self.gs.build("V", zen, azi, vec![mag as f32]);

        self.alt_az_pa = self.probe.alt_az_parallactic(&self.obs);
        self.gs.rotate_rays(self.alt_az_pa.2);

        self.gs.through(&mut self.gmt).xpupil();
        self.sensor.calibrate(&mut self.gs, 0.9);
        self
    }
    pub fn propagate_src(&mut self) {
        self.gs
            .through(&mut self.gmt)
            .xpupil()
            .through(&mut self.sensor);
    }
    pub fn update(&mut self, inc_secs: f64, gstate: &ceo::GmtState) -> &mut Self {
        self.obs.add_seconds(inc_secs);
        self.gmt.update(gstate);

        let alt_az_pa = self.probe.alt_az_parallactic(&self.obs);
        self.gs.rotate_rays(alt_az_pa.2);

        self.gs
            .through(&mut self.gmt)
            .xpupil()
            .through(&mut self.sensor);
        self
    }
    pub fn local(&mut self) -> (f64, f64, f64) {
        let (z, a) = self
            .probe
            .local_polar_arcmin_deg(&self.telescope, &self.obs);
        let alt_az_pa = self.probe.alt_az_parallactic(&self.obs);
        (z, a, (self.alt_az_pa.2 - alt_az_pa.2) / DEG2RAD)
    }
    pub fn calibrate(&mut self, progress: Option<ProgressBar>) -> Array2<f32> {
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

pub struct OpticalPathToDSH48 {
    pub obs: Observation,
    pub telescope: SkyCoordinates,
    pub probe: SkyCoordinates,
    pub gmt: ceo::Gmt,
    pub gs: ceo::Source,
    pub sensor: ceo::ShackHartmann,
    pub atm: ceo::Atmosphere,
    pub probe_id: i32,
    tel_alt_az_pa: (f64, f64, f64),
    za: (f64, f64),
}
impl OpticalPathToDSH48 {
    pub fn new(
        datetime: &str,
        telescope_radec_deg: (f64, f64),
        probe_radec_deg: (f64, f64),
        probe_id: i32,
    ) -> OpticalPathToDSH48 {
        let n_side_lenslet = 48;
        let n_px_lenslet = 16;
        //let n_px = n_side_lenslet*16 + 1;
        let pupil_size = 25.5;
        let m1_n_mode = 27;
        let d = pupil_size / n_side_lenslet as f64;
        OpticalPathToDSH48 {
            obs: Observation::new(datetime, GMT_LAT, GMT_LONG),
            telescope: SkyCoordinates::new(telescope_radec_deg.0, telescope_radec_deg.1),
            probe: SkyCoordinates::new(probe_radec_deg.0, probe_radec_deg.1),
            gmt: ceo::Gmt::new(m1_n_mode, None),
            gs: ceo::Source::empty(),
            sensor: ceo::ShackHartmann::new(1, n_side_lenslet, n_px_lenslet, d),
            atm: ceo::Atmosphere::new(),
            probe_id,
            tel_alt_az_pa: (0., 0., 0.),
            za: (0., 0.),
        }
    }
    pub fn build(
        &mut self,
        n_px_framelet: i32,
        n_px_imagelet: Option<i32>,
        osf: Option<i32>,
        mag: f64,
    ) -> &mut Self {
        self.gmt.build();
        self.sensor.build(n_px_framelet, n_px_imagelet, osf);
        self.gs = self.sensor.new_guide_stars();
        let q = self.probe.local(&self.telescope, &self.obs);
        let z = q[0].hypot(q[1]) as f32;
        let a = q[1].atan2(q[0]) as f32;
        let zen: Vec<f32> = vec![z]; //6.0*f32::consts::PI/180./60.];
        let azi: Vec<f32> = vec![a]; //2. * (self.probe_id as f32) * f32::consts::PI / 3.];
        self.gs.build("V", zen, azi, vec![mag as f32]);

        self.tel_alt_az_pa = self.telescope.alt_az_parallactic(&self.obs);
        self.za = self
            .probe
            .local_polar_arcmin_deg(&self.telescope, &self.obs);
        self.gs.rotate_rays(self.tel_alt_az_pa.2);
        self.gs.set_fwhm(3.16);

        self.gs.through(&mut self.gmt).xpupil();
        self.sensor.calibrate(&mut self.gs, 0.9);
        self
    }
    pub fn build_atmosphere(&mut self, fullpath_to_phasescreens: &str) {
        self.atm.load_from_json(fullpath_to_phasescreens).unwrap();
    }
    pub fn propagate_src(&mut self) {
        self.gs
            .through(&mut self.gmt)
            .xpupil()
            .through(&mut self.sensor);
    }
    pub fn update(&mut self, inc_secs: f64, gstate: &ceo::GmtState) -> &mut Self {
        self.obs.add_seconds(inc_secs);
        self.gmt.update(gstate);

        let alt_az_pa = self.probe.alt_az_parallactic(&self.obs);
        /*
        self.gs.rotate_rays(alt_az_pa.2);

        let atm_sampling = 1e-2;
        let n = (inc_secs / atm_sampling) as u32;
        for k in 0..n {
            self.gs
                .through(&mut self.gmt)
                .xpupil()
                .through(&mut self.atm)
                .through(&mut self.sensor);
            self.atm.secs += atm_sampling;
        }
        */
        self.gs
            .through(&mut self.gmt)
            .xpupil()
            .through(&mut self.sensor);
        self
    }
    pub fn local(&mut self) -> (f64, f64, f64, f64) {
        let (z, a) = self
            .probe
            .local_polar_arcmin_deg(&self.telescope, &self.obs);
        let alt_az_pa = self.telescope.alt_az_parallactic(&self.obs);
        (
            z,
            a,
            self.za.1 - a,
            (self.tel_alt_az_pa.2 - alt_az_pa.2) / DEG2RAD,
        )
    }
}
impl Drop for OpticalPathToDSH48 {
    fn drop(&mut self) {
        drop(&mut self.gmt);
        drop(&mut self.gs);
        drop(&mut self.sensor);
    }
}
