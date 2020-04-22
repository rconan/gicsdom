use std::f64;

/// GMT observatory latitude
pub const GMT_LAT: f64 = -29.049;
/// GMT observatory longitude
pub const GMT_LONG: f64 = -70.682;

#[derive(Clone, Debug)]
/// A module to manage time system and conversion
///
/// [Practical Astronomy, P. Duffet-Smith and J. Zwart][1]
///
/// [1]: https://archive.org/details/Practical_Astronomy_with_your_Calculator_or_Spreadsheet_4th_edition_by_Peter_Du
pub struct Time {
    /// The year
    pub y: i32,
    /// The month
    pub m: u8,
    /// The day
    pub d: u8,
    /// The hour
    pub h: u8,
    /// The minutes
    pub mn: u8,
    /// The seconds
    pub s: u8,
    /// The fractional second
    pub fs: f64,
}
impl Time {
    /// Converts a date and a time in UTC to a `Time`
    pub fn from_date_utc(y: i32, m: u8, d: u8, h: u8, mn: u8, s: u8) -> Self {
        Time {
            y,
            m,
            d,
            h,
            mn,
            s,
            fs: 0.0,
        }
    }
    /// Converts a date and a time in UTC down to the fractional second to a `Time`
    pub fn from_date_utc_fs(y: i32, m: u8, d: u8, h: u8, mn: u8, s: u8, fs: f64) -> Self {
        Time {
            y,
            m,
            d,
            h,
            mn,
            s,
            fs,
        }
    }
    /// Converts a time in UTC to a `Time`
    pub fn from_utc(h: u8, mn: u8, s: u8) -> Self {
        Time {
            y: 0,
            m: 0,
            d: 0,
            h,
            mn,
            s: s,
            fs: 0.0,
        }
    }
    /// Returns the date and time as a `String` following ISO 8601 standard
    pub fn datetime(&self) -> String {
        format!(
            "{:4}-{:02}-{:02}T{:02}:{:02}:{:06.3}",
            self.y,
            self.m,
            self.d,
            self.h,
            self.mn,
            self.s as f64 + self.fs
        )
    }
    /// Returns the decimal days
    pub fn decimal_days(&self) -> f64 {
        self.d as f64
            + (self.h as f64 + self.mn as f64 * 60. + (self.s as f64 + self.fs as f64) * 3600.)
                / 24.
    }
    /// Returns the decimal hours
    pub fn decimal_hours(&self) -> f64 {
        (self.h as f64) + ((self.mn as f64) + (self.s as f64 + self.fs as f64) / 60.0) / 60.0
    }
    /// Returns the Julian date
    pub fn julian_date(&self) -> f64 {
        if self.y < 1583 {
            unimplemented!("Julian date is implemented for years after 1582!");
        }
        let (yp, mp) = if self.m == 1 || self.m == 2 {
            (self.y as f64 - 1., self.m as f64 + 12.)
        } else {
            (self.y as f64, self.m as f64)
        };
        let a = (yp / 100.).trunc();
        let b = 2. - a + (a / 4.).trunc();
        let c = if yp < 0.0 {
            (365.25 * yp - 0.75).trunc()
        } else {
            (365.25 * yp).trunc()
        };
        let d = (30.6001 * (mp + 1.)).trunc();
        b + c + d + self.decimal_days() + 1720994.5
    }
    /// Returns the Julian date at midnight
    pub fn julian_date_at_midnight(&self) -> f64 {
        if self.y < 1583 {
            unimplemented!("Julian date is implemented for years after 1582!");
        }
        let (yp, mp) = if self.m == 1 || self.m == 2 {
            (self.y as f64 - 1., self.m as f64 + 12.)
        } else {
            (self.y as f64, self.m as f64)
        };
        let a = (yp / 100.).trunc();
        let b = 2. - a + (a / 4.).trunc();
        let c = if yp < 0.0 {
            (365.25 * yp - 0.75).trunc()
        } else {
            (365.25 * yp).trunc()
        };
        let d = (30.6001 * (mp + 1.)).trunc();
        b + c + d + self.d as f64 + 1720994.5
    }
    /// Returns the Greenwich sidereal time
    pub fn greenwich_sidereal_time(&self) -> f64 {
        let s = self.julian_date_at_midnight() - 2451545.0;
        let t = s / 36525.;
        let t0 = (6.697374558 + (2400.051336 + 0.000025862 * t) * t).rem_euclid(24.0);
        let ut = self.decimal_hours();
        (ut * 1.002737909 + t0).rem_euclid(24.0)
    }
    /// Returns the local sidereal time
    pub fn local_sidereal_time(&self, longitude_deg: f64) -> f64 {
        self.greenwich_sidereal_time() + longitude_deg / 15.
    }
    /// Adds seconds to `Time`
    pub fn add_seconds(&mut self, secs: f64) {
        let mut s = self.decimal_hours() * 3600. + secs;
        self.h = (s / 3600.).trunc() as u8;
        if self.h >= 24 {
            let d = (self.h / 24) as u8;
            self.d += d;
            self.h -= (self.d * 24) as u8;
        }
        s -= self.h as f64 * 3600.;
        self.mn = (s / 60.).trunc() as u8;
        s -= self.mn as f64 * 60.;
        self.s = s.round() as u8;
        self.fs = s - self.s as f64;
    }
}
#[derive(Clone, Debug)]
/// Specifications of an observation
pub struct Observation {
    /// site latitude [radian]
    pub latitude: f64,
    /// site longitude [radian]
    pub longitude: f64,
    /// observation date and time
    pub utc: Time,
    /// `SkyCoordinates` of the observation celestial target
    pub object: SkyCoordinates,
    /// observation sampling time [s]
    pub sampling_time: f64,
    /// observation duration [s]
    pub duration: f64,
    /// observation number of step since the beginning of the observation; step=[0,duration/sampling_time-1]
    pub step: u32,
    /// a flag raised when the observation is complete i.e. step=duration/sampling_time-1
    pub ended: bool,
}
impl Observation {
    /// Returns an observation roster
    ///
    /// # Arguments
    ///
    /// * `latitude_deg` - the observatory latitude [degree]
    /// * `longitude_deg` - the observatory longitude [degree]
    /// * `utc` - the observation `Time`
    /// * `object` - celestial target `SkyCoordinates`
    /// * `sampling_time` - observation sampling time [s]
    /// * `duration` - observation duration [s]
    pub fn new(
        latitude_deg: f64,
        longitude_deg: f64,
        utc: Time,
        object: SkyCoordinates,
        sampling_time: f64,
        duration: f64,
    ) -> Observation {
        Observation {
            latitude: latitude_deg.to_radians(),
            longitude: longitude_deg.to_radians(),
            utc,
            object,
            sampling_time,
            duration,
            step: 0,
            ended: false,
        }
    }
    pub fn from_date_utc(
        latitude_deg: f64,
        longitude_deg: f64,
        utc: Time,
        object: SkyCoordinates,
        sampling_time: f64,
        duration: f64,
    ) -> Observation {
        Observation {
            latitude: latitude_deg.to_radians(),
            longitude: longitude_deg.to_radians(),
            utc,
            object,
            sampling_time,
            duration,
            step: 0,
            ended: false,
        }
    }
    /// Returns the local sidereal time
    pub fn lst(&self) -> f64 {
        self.utc.local_sidereal_time(self.longitude.to_degrees())
    }
}
impl Iterator for Observation {
    type Item = u32;
    /// Increments `step` and adds `sampling_time` seconds to `utc`
    fn next(&mut self) -> Option<Self::Item> {
        self.step += 1;
        let s = self.step as f64 * self.sampling_time;
        self.utc.add_seconds(self.sampling_time);
        if s <= self.duration {
            Some(self.step)
        } else {
            self.ended = true;
            None
        }
    }
}

#[derive(Clone, Debug)]
/// Celestial system of coordinates and system conversion
pub struct SkyCoordinates {
    /// right ascension and declination [degree]
    pub radec: (f64, f64),
}
/// Transforms Cartesian coordinates `vec![x,y]` into polar coordinates `vec![r,o]`
pub fn pol2cart(v: Vec<f64>) -> Vec<f64> {
    let v_pol: Vec<f64> = vec![v[0].hypot(v[1]), v[1].atan2(v[0])];
    v_pol
}
impl SkyCoordinates {
    /// Returns a new sky coordinate system using equatorial coordinates right ascension [degree] and declination [degree]
    pub fn new(ra_deg: f64, dec_deg: f64) -> SkyCoordinates {
        SkyCoordinates {
            radec: (ra_deg.to_radians(), dec_deg.to_radians()),
        }
    }
    /// Returns equatorial coordinates right ascension [degree] and declination [degree]
    pub fn radec_deg(&self) -> (f64, f64) {
        (self.radec.0.to_degrees(), self.radec.1.to_degrees())
    }
    /// Returns the hour angle for a given `Observation`
    pub fn hour_angle(&self, obs: &Observation) -> f64 {
        let ha = obs.lst() - self.radec.0.to_degrees() / 15.;
        if ha < 0. {
            ha + 24.
        } else {
            ha
        }
    }
    /// Returns the horizon coordinates altitude [rd], azimuth [rd] and parallactic angle [rd] for a given `Observation`
    pub fn alt_az_parallactic(&self, obs: &Observation) -> (f64, f64, f64) {
        let sin_cos_ha = (15. * self.hour_angle(obs)).to_radians().sin_cos();
        let sin_alt = self.radec.1.sin() * obs.latitude.sin()
            + self.radec.1.cos() * obs.latitude.cos() * sin_cos_ha.1;
        let alt = sin_alt.asin();
        let cos_az = obs.latitude.sin() * sin_cos_ha.1 - obs.latitude.cos() * self.radec.1.tan();
        //        let cos_az = sin_alt * obs.latitude.sin() - self.radec.1.sin();
        let az = sin_cos_ha.0.atan2(cos_az);
        //        let az = cos_az.atan2(sin_cos_ha.0);
        let cos_p = obs.latitude.tan() * self.radec.1.cos() - self.radec.1.sin() * sin_cos_ha.1;
        let p = sin_cos_ha.0.atan2(cos_p);
        (alt, az, p)
    }
    /// Returns the horizon coordinates altitude [rd] and azimuth [rd]
    pub fn altaz(&self, obs: &Observation) -> (f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0, aap.1)
    }
    /// Returns the horizon coordinates altitude [rd] and azimuth [degree]
    pub fn altaz_deg(&self, obs: &Observation) -> (f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0.to_degrees(), aap.1.to_degrees())
    }
    /// Returns the horizon coordinates altitude [rd]]
    pub fn alt(&self, obs: &Observation) -> f64 {
        self.altaz(obs).0
    }
    /// Returns the zenith angle [rd]]
    pub fn za(&self, obs: &Observation) -> f64 {
        f64::consts::FRAC_PI_2 - self.alt(obs)
    }
    /// Returns the horizon coordinates altitude [degree], azimuth [degree] and parallactic angle [degree] for a given `Observation`
    pub fn alt_az_parallactic_deg(&self, obs: &Observation) -> (f64, f64, f64) {
        let aap = self.alt_az_parallactic(obs);
        (aap.0.to_degrees(), aap.1.to_degrees(), aap.2.to_degrees())
    }
    /// Returns the local coordinates for a given `Observation`.
    /// The local coordinates are given be the transformation of the horizon coordinates into the horizon coordinates of the object of the `Observation`
    pub fn local(&self, obs: &Observation) -> Vec<f64> {
        let altaz = self.altaz(obs);
        let sc_alt = altaz.0.sin_cos();
        let sc_az = altaz.1.sin_cos();
        let (x, y, z) = (sc_az.1 * sc_alt.1, sc_az.0 * sc_alt.1, sc_alt.0);
        let ref_altaz = obs.object.altaz(obs);
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
    /// Returns the local polar coordinates for a given `Observation`.
    pub fn local_polar(&self, obs: &Observation) -> (f64, f64) {
        let u = self.local(obs);
        let v = pol2cart(u);
        (v[0], v[1])
    }
    /// Returns the local polar coordinates in arcmin for the radius and in degree for the azimuth, for a given `Observation`.
    pub fn local_polar_arcmin_deg(&self, obs: &Observation) -> (f64, f64) {
        let (r, o) = self.local_polar(obs);
        (60. * r.to_degrees(), o.to_degrees())
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_decimal_days() {
        let t = Time::from_date_utc(2009, 6, 19, 18, 0, 0);
        assert_eq!(t.decimal_days(), 19.75)
    }

    #[test]
    fn time_decimal_hours() {
        let t = Time::from_date_utc_fs(1980, 4, 22, 14, 36, 51, 0.67);
        assert!((t.decimal_hours() - 14.614353).abs() < 1e-6)
    }

    #[test]
    fn time_julian_date() {
        let t = Time::from_date_utc(2009, 6, 19, 18, 0, 0);
        assert_eq!(t.julian_date(), 2455002.25)
    }

    #[test]
    fn time_gst() {
        let t = Time::from_date_utc_fs(1980, 4, 22, 14, 36, 51, 0.67);
        assert!((t.greenwich_sidereal_time() - 4.668120).abs() < 1e-6)
    }

    #[test]
    fn time_lst() {
        let t = Time::from_date_utc_fs(1980, 4, 22, 14, 36, 51, 0.67);
        assert!((t.local_sidereal_time(-64.) - 0.401453).abs() < 1e-6)
    }

    #[test]
    fn sky_hour_angle() {
        let t = Time::from_date_utc_fs(1980, 4, 22, 18, 36, 51, 0.67);
        let ra = Time::from_utc(18, 32, 21);
        let s = SkyCoordinates::new(ra.decimal_hours() * 15., 0.0);
        let obs = Observation::from_date_utc(0., -64., t, s, 1., 1.);
        assert!((obs.object.hour_angle(&obs) - 9.873237).abs() < 1e-6)
    }

    #[test]
    fn sky_altaz() {
        let t = Time::from_date_utc(2012, 06, 10, 4, 1, 3);
        let s = SkyCoordinates::new(266.62156258190726, -27.776114821065107);
        let obs = Observation::from_date_utc(GMT_LAT, GMT_LONG, t, s, 1., 1.);
        let (alt, az, p) = obs.object.alt_az_parallactic(&obs);
        println!(
            "alt={},az={},p={}",
            alt.to_degrees(),
            az.to_degrees(),
            p.to_degrees()
        );
    }
}
