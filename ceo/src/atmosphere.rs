use serde::Deserialize;
use std::ffi::CString;
use std::{f32, mem};

use super::ceo_bindings::atmosphere;
use super::{element, Propagation, Source, CEO};

#[derive(Deserialize, Debug)]
struct GmtAtmosphere {
    r0: f32,
    #[serde(rename = "L0")]
    l_not: f32,
    #[serde(rename = "L")]
    length: f32,
    #[serde(rename = "lower_case")]
    nxy_pupil: i32,
    fov: f32,
    duration: f32,
    #[serde(rename = "lower_case")]
    n_duration: i32,
    filename: String,
    #[serde(rename = "lower_case")]
    seed: i32,
}

pub struct Atmosphere {
    _c_: atmosphere,
    pub r0_at_zenith: f64,
    pub oscale: f64,
    pub zenith_angle: f64,
    pub secs: f64,
    //filename: String,
    //k_duration: i32,
    propagate_ptr: fn(&mut Atmosphere, &mut Source, f32),
}
/// ## `Atmosphere` builder
impl CEO<element::ATMOSPHERE> {
    /// Create a new `Atmosphere` builder
    pub fn new() -> CEO<element::ATMOSPHERE> {
        CEO {
            args: element::ATMOSPHERE::default(),
        }
    }
    /// Set r0 value taken at pointing the zenith in meters
    pub fn set_r0_at_zenith(mut self, r0_at_zenith: f64) -> Self {
        self.args.r0_at_zenith = r0_at_zenith;
        self
    }
    /// Set outer scale value in meters
    pub fn set_oscale(mut self, oscale: f64) -> Self {
        self.args.oscale = oscale;
        self
    }
    /// Set zenith angle value in radians
    pub fn set_zenith_angle(mut self, zenith_angle: f64) -> Self {
        self.args.zenith_angle = zenith_angle;
        self
    }
    pub fn set_turbulence_profile(mut self, turbulence: element::TurbulenceProfile) -> Self {
        self.args.turbulence = turbulence;
        self
    }
    /// Set parameters for atmosphere ray tracing
    pub fn set_ray_tracing(
        mut self,
        width: f32,
        n_width_px: i32,
        field_size: f32,
        duration: f32,
        filepath: Option<String>,
        n_duration: Option<i32>,
    ) -> Self {
        self.args.ray_tracing = Some(element::RayTracing {
            width,
            n_width_px,
            field_size,
            duration,
            filepath,
            n_duration,
        });
        self
    }
    /// Build the `Atmosphere`
    pub fn build(mut self) -> Atmosphere {
        let mut atm = Atmosphere {
            _c_: unsafe { mem::zeroed() },
            r0_at_zenith: self.args.r0_at_zenith,
            oscale: self.args.oscale,
            zenith_angle: self.args.zenith_angle,
            secs: 0.0,
            //filename: String::new(),
            //k_duration: 0,
            propagate_ptr: |_, _, _| (),
        };
        let r0 = (atm.r0_at_zenith.powf(-5.0 / 3.0) / atm.zenith_angle.cos()).powf(-3.0 / 5.0);
        match self.args.ray_tracing {
            None => unsafe {
                atm._c_.setup(
                    r0 as f32,
                    self.args.oscale as f32,
                    self.args.turbulence.n_layer as i32,
                    self.args.turbulence.altitude.as_mut_ptr(),
                    self.args.turbulence.xi0.as_mut_ptr(),
                    self.args.turbulence.wind_speed.as_mut_ptr(),
                    self.args.turbulence.wind_direction.as_mut_ptr(),
                );
                atm.propagate_ptr = |a, s, t| {
                    let n_xy = s.pupil_sampling;
                    let d_xy = (s.pupil_size / (n_xy - 1) as f64) as f32;
                    a._c_
                        .get_phase_screen4(s.as_raw_mut_ptr(), d_xy, n_xy, d_xy, n_xy, t);
                };
            },
            Some(rtc) => match rtc.filepath {
                Some(file) => unsafe {
                    atm._c_.setup2(
                        r0 as f32,
                        self.args.oscale as f32,
                        self.args.turbulence.n_layer as i32,
                        self.args.turbulence.altitude.as_mut_ptr(),
                        self.args.turbulence.xi0.as_mut_ptr(),
                        self.args.turbulence.wind_speed.as_mut_ptr(),
                        self.args.turbulence.wind_direction.as_mut_ptr(),
                        rtc.width,
                        rtc.n_width_px,
                        rtc.field_size,
                        rtc.duration,
                        CString::new(file.to_owned().into_bytes())
                            .unwrap()
                            .into_raw(),
                        rtc.n_duration.unwrap_or(1),
                    );
                    atm.propagate_ptr = |a, s, t| {
                        let n_xy = s.pupil_sampling;
                        let d_xy = (s.pupil_size / (n_xy - 1) as f64) as f32;
                        a._c_
                            .rayTracing1(s.as_raw_mut_ptr(), d_xy, n_xy, d_xy, n_xy, t);
                    };
                },
                None => unsafe {
                    atm._c_.setup1(
                        r0 as f32,
                        self.args.oscale as f32,
                        self.args.turbulence.n_layer as i32,
                        self.args.turbulence.altitude.as_mut_ptr(),
                        self.args.turbulence.xi0.as_mut_ptr(),
                        self.args.turbulence.wind_speed.as_mut_ptr(),
                        self.args.turbulence.wind_direction.as_mut_ptr(),
                        rtc.width,
                        rtc.n_width_px,
                        rtc.field_size,
                        rtc.duration,
                    );
                    atm.propagate_ptr = |a, s, t| {
                        let n_xy = s.pupil_sampling;
                        let d_xy = (s.pupil_size / (n_xy - 1) as f64) as f32;
                        a._c_
                            .rayTracing1(s.as_raw_mut_ptr(), d_xy, n_xy, d_xy, n_xy, t);
                    };
                },
            },
        }
        atm
    }
}
impl Atmosphere {
    pub fn new() -> Atmosphere {
        Atmosphere {
            _c_: unsafe { mem::zeroed() },
            r0_at_zenith: 0.16,
            oscale: 25.5,
            zenith_angle: 0.0,
            secs: 0.0,
            //filename: String::new(),
            //k_duration: 0,
            propagate_ptr: |_, _, _| (),
        }
    }
    pub fn build(
        &mut self,
        r_not: f32,
        l_not: f32,
        n_layer: i32,
        mut altitude: Vec<f32>,
        mut xi0: Vec<f32>,
        mut wind_speed: Vec<f32>,
        mut wind_direction: Vec<f32>,
    ) -> &mut Self {
        unsafe {
            self._c_.setup(
                r_not,
                l_not,
                n_layer,
                altitude.as_mut_ptr(),
                xi0.as_mut_ptr(),
                wind_speed.as_mut_ptr(),
                wind_direction.as_mut_ptr(),
            );
        }
        self.propagate_ptr = |a, s, t| unsafe {
            let n_xy = s.pupil_sampling;
            let d_xy = (s.pupil_size / (n_xy - 1) as f64) as f32;
            a._c_
                .get_phase_screen4(s.as_raw_mut_ptr(), d_xy, n_xy, d_xy, n_xy, t);
        };
        self
    }
    pub fn as_raw_mut_ptr(&mut self) -> &mut atmosphere {
        &mut self._c_
    }
    pub fn raytrace_build(
        &mut self,
        r_not: f32,
        l_not: f32,
        n_layer: i32,
        mut altitude: Vec<f32>,
        mut xi0: Vec<f32>,
        mut wind_speed: Vec<f32>,
        mut wind_direction: Vec<f32>,
        width: f32,
        n_width_px: i32,
        field_size: f32,
        duration: f32,
        filepath: Option<&str>,
        n_duration: Option<i32>,
    ) -> &mut Self {
        match filepath {
            Some(file) => unsafe {
                self._c_.setup2(
                    r_not,
                    l_not,
                    n_layer,
                    altitude.as_mut_ptr(),
                    xi0.as_mut_ptr(),
                    wind_speed.as_mut_ptr(),
                    wind_direction.as_mut_ptr(),
                    width,
                    n_width_px,
                    field_size,
                    duration,
                    CString::new(file.to_owned().into_bytes())
                        .unwrap()
                        .into_raw(),
                    n_duration.unwrap_or(1),
                );
            },
            None => unsafe {
                self._c_.setup1(
                    r_not,
                    l_not,
                    n_layer,
                    altitude.as_mut_ptr(),
                    xi0.as_mut_ptr(),
                    wind_speed.as_mut_ptr(),
                    wind_direction.as_mut_ptr(),
                    width,
                    n_width_px,
                    field_size,
                    duration,
                );
            },
        }
        self.propagate_ptr = |a, s, t| unsafe {
            let n_xy = s.pupil_sampling;
            let d_xy = (s.pupil_size / (n_xy - 1) as f64) as f32;
            a._c_
                .rayTracing1(s.as_raw_mut_ptr(), d_xy, n_xy, d_xy, n_xy, t);
        };
        self
    }
    pub fn gmt_build(&mut self, r_not: f32, l_not: f32) -> &mut Self {
        unsafe {
            self._c_.gmt_setup4(r_not, l_not, 2020);
        }
        self
    }
    pub fn set_r0(&mut self, new_r0: f64) {
        self._c_.r0 = new_r0 as f32;
    }
    pub fn reset(&mut self) {
        unsafe {
            self._c_.reset();
        }
    }
}
impl Drop for Atmosphere {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Atmosphere {
    fn time_propagate(&mut self, secs: f64, src: &mut Source) -> &mut Self {
        (self.propagate_ptr)(self, src, secs as f32);
        self
    }
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        self.time_propagate(self.secs, src)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atmosphere_new() {
        crate::ceo!(element::ATMOSPHERE);
    }
}
