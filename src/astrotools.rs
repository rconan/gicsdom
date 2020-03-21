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
        lst.to_radians()
    }
    pub fn local_sidereal_time_deg(&self) -> f64 {
        self.local_sidereal_time().to_degrees()
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
