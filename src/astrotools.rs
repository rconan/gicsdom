use std::f64;
use hifitime::Epoch;

pub const GMT_LAT: f64 = -29.049;
pub const GMT_LONG: f64 = -70.682;

#[derive(Clone, Debug)]
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
            latitude: latitude_deg.to_radians(),
            longitude: longitude_deg.to_radians(),
            epoch: Epoch::from_gregorian_utc_str(_datetime).unwrap(),
        }
    }
    pub fn latitude_deg(&self) -> f64 {
        self.latitude.to_degrees()
    }
    pub fn longitude_deg(&self) -> f64 {
        self.longitude.to_degrees()
    }
    pub fn julian_day(&self) -> f64 {
        self.epoch.as_jde_utc_days()
    }
    pub fn j2000_day(&self) -> f64 {
        self.julian_day() - 2451545.0
    }
    pub fn greenwich_sideral_time(&self) -> f64 {
        let (y, m, d, h, min, sec) = self.epoch.as_gregorian_utc();
        let jd = Epoch::from_gregorian_utc_at_midnight(y,m,d).as_jde_utc_days(); 
        let s = jd - 2451545.0;
        let t = s/36525.;
        let t0 = (6.697374558+(2400.051336+0.000025862*t)*t).rem_euclid(24.0);
        let ut =(h as f64) + ((min as f64) + (sec as f64) / 60.0) / 60.0;
        (ut*1.002737909 + t0).rem_euclid(24.0)
    }
    pub fn local_sidereal_time(&self) -> f64 {
        (self.greenwich_sideral_time()*15.+self.longitude.to_degrees()).to_radians()
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
            radec: (ra_deg.to_radians(), dec_deg.to_radians()),
        }
    }
    pub fn radec_deg(&self) -> (f64, f64) {
        (self.radec.0.to_degrees(), self.radec.1.to_degrees())
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
        (aap.0.to_degrees(), aap.1.to_degrees())
    }
    pub fn alt(&self, obs: &Observation) -> f64 {
        self.altaz(obs).0
    }
    pub fn za(&self, obs: &Observation) -> f64 {
        f64::consts::FRAC_PI_2 - self.alt(obs)
    }
    pub fn alt_az_parallactic_deg(&self, obs: &Observation) -> (f64, f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0.to_degrees(), aap.1.to_degrees(), aap.2.to_degrees())
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
    pub fn local_polar(&self, reference: &SkyCoordinates, obs: &Observation) -> (f64, f64) {
        let u = self.local(reference, obs);
        let v = pol2cart(u);
        (v[0], v[1])
    }
    pub fn local_polar_arcmin_deg(
        &self,
        reference: &SkyCoordinates,
        obs: &Observation,
    ) -> (f64, f64) {
        let (r, o) = self.local_polar(reference, obs);
        (60.*r.to_degrees(), o.to_degrees())
    }
}
