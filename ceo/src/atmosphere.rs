use log;
use serde::Deserialize;
use std::ffi::CString;
use std::{f32, mem};

use super::ceo_bindings::atmosphere;
use super::{Builder, Propagation, Source};

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

#[derive(Debug, Clone)]
#[doc(hidden)]
pub struct TurbulenceProfile {
    pub n_layer: usize,
    pub altitude: Vec<f32>,
    pub xi0: Vec<f32>,
    pub wind_speed: Vec<f32>,
    pub wind_direction: Vec<f32>,
}
impl Default for TurbulenceProfile {
    fn default() -> Self {
        TurbulenceProfile {
            n_layer: 7,
            altitude: [25.0, 275.0, 425.0, 1_250.0, 4_000.0, 8_000.0, 13_000.0].to_vec(),
            xi0: [0.1257, 0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751].to_vec(),
            wind_speed: [5.6540, 5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187].to_vec(),
            wind_direction: [0.0136, 0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462].to_vec(),
        }
    }
}
#[derive(Debug, Clone)]
#[doc(hidden)]
pub struct RayTracing {
    pub width: f32,
    pub n_width_px: i32,
    pub field_size: f32,
    pub duration: f32,
    pub filepath: Option<String>,
    pub n_duration: Option<i32>,
}
/// [`CEO`](../struct.CEO.html#impl-6) [`Atmosphere`](../struct.Atmosphere.html) builder type
#[derive(Debug, Clone)]
pub struct ATMOSPHERE {
    pub r0_at_zenith: f64,
    pub oscale: f64,
    pub zenith_angle: f64,
    pub turbulence: TurbulenceProfile,
    pub ray_tracing: Option<RayTracing>,
}
/// Default properties:
///  * r0           : 16cm
///  * L0           : 25m
///  * zenith angle : 30 degrees
///  * turbulence profile:
///    * n_layer        : 7
///    * altitude       : [25.0, 275.0, 425.0, 1250.0, 4000.0, 8000.0, 13000.0] m
///    * xi0            : [0.1257, 0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751]
///    * wind speed     : [5.6540, 5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187] m/s
///    * wind direction : [0.0136, 0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462] rd
/// * ray tracing : none
impl Default for ATMOSPHERE {
    fn default() -> Self {
        ATMOSPHERE {
            r0_at_zenith: 0.16,
            oscale: 25.,
            zenith_angle: 30_f64.to_radians(),
            turbulence: TurbulenceProfile::default(),
            ray_tracing: None,
        }
    }
}
/// ## `Atmosphere` builder
impl ATMOSPHERE {
    /// Set r0 value taken at pointing the zenith in meters
    pub fn set_r0_at_zenith(self, r0_at_zenith: f64) -> Self {
        Self {
            r0_at_zenith,
            ..self
        }
    }
    /// Set outer scale value in meters
    pub fn set_oscale(self, oscale: f64) -> Self {
        Self { oscale, ..self }
    }
    /// Set zenith angle value in radians
    pub fn set_zenith_angle(self, zenith_angle: f64) -> Self {
        Self {
            zenith_angle,
            ..self
        }
    }
    /// Set the turbulence profile
    pub fn set_turbulence_profile(self, turbulence: TurbulenceProfile) -> Self {
        Self { turbulence, ..self }
    }
    /// Set a single turbulence layer
    pub fn set_single_turbulence_layer(
        self,
        altitude: f32,
        wind_speed: Option<f32>,
        wind_direction: Option<f32>,
    ) -> Self {
        Self {
            turbulence: TurbulenceProfile {
                n_layer: 1,
                altitude: vec![altitude],
                xi0: vec![1f32],
                wind_speed: vec![wind_speed.unwrap_or(0f32)],
                wind_direction: vec![wind_direction.unwrap_or(0f32)],
            },
            ..self
        }
    }
    /// Remove a turbulence layer specifield by its zero based index
    pub fn remove_turbulence_layer(self, layer_idx: usize) -> Self {
        let mut turbulence = self.turbulence;
        turbulence.n_layer -= 1;
        turbulence.altitude.remove(layer_idx);
        turbulence.xi0.remove(layer_idx);
        turbulence.wind_speed.remove(layer_idx);
        turbulence.wind_direction.remove(layer_idx);
        Self { turbulence, ..self }
    }
    /// Set parameters for atmosphere ray tracing
    pub fn set_ray_tracing(
        self,
        width: f32,
        n_width_px: i32,
        field_size: f32,
        duration: f32,
        filepath: Option<String>,
        n_duration: Option<i32>,
    ) -> Self {
        Self {
            ray_tracing: Some(RayTracing {
                width,
                n_width_px,
                field_size,
                duration,
                filepath,
                n_duration,
            }),
            ..self
        }
    }
}
impl Builder for ATMOSPHERE {
    type Component = Atmosphere;
    /// Build the `Atmosphere`
    fn build(self) -> Atmosphere {
        let mut atm = Atmosphere {
            _c_: unsafe { mem::zeroed() },
            r0_at_zenith: self.r0_at_zenith,
            oscale: self.oscale,
            zenith_angle: self.zenith_angle,
            secs: 0.0,
            //filename: String::new(),
            //k_duration: 0,
            propagate_ptr: |_, _, _| (),
        };
        let secz = 1f64 / atm.zenith_angle.cos();
        let r0 = (atm.r0_at_zenith.powf(-5.0 / 3.0) * secz).powf(-3.0 / 5.0);
        log::info!(
            "Atmosphere r0 at {:.1}degree from zenith: {:.3}m",
            atm.zenith_angle.to_degrees(),
            r0
        );
        let mut altitude = self
            .turbulence
            .altitude
            .iter()
            .map(|x| *x as f32 * secz as f32)
            .collect::<Vec<f32>>();
        let mut wind_speed = self
            .turbulence
            .wind_speed
            .iter()
            .map(|x| *x as f32 / secz as f32)
            .collect::<Vec<f32>>();
        let mut xi0 = self.turbulence.xi0;
        let mut wind_direction = self.turbulence.wind_direction;
        match &self.ray_tracing {
            None => unsafe {
                atm._c_.setup(
                    r0 as f32,
                    self.oscale as f32,
                    self.turbulence.n_layer as i32,
                    altitude.as_mut_ptr(),
                    xi0.as_mut_ptr(),
                    wind_speed.as_mut_ptr(),
                    wind_direction.as_mut_ptr(),
                );
                atm.propagate_ptr = |a, s, t| {
                    let n_xy = s.pupil_sampling;
                    let d_xy = (s.pupil_size / (n_xy - 1) as f64) as f32;
                    a._c_
                        .get_phase_screen4(s.as_raw_mut_ptr(), d_xy, n_xy, d_xy, n_xy, t);
                };
            },
            Some(rtc) => match &rtc.filepath {
                Some(file) => unsafe {
                    log::info!("Looking up phase screen from file {}", file);
                    atm._c_.setup2(
                        r0 as f32,
                        self.oscale as f32,
                        self.turbulence.n_layer as i32,
                        altitude.as_mut_ptr(),
                        xi0.as_mut_ptr(),
                        wind_speed.as_mut_ptr(),
                        wind_direction.as_mut_ptr(),
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
                        self.oscale as f32,
                        self.turbulence.n_layer as i32,
                        altitude.as_mut_ptr(),
                        xi0.as_mut_ptr(),
                        wind_speed.as_mut_ptr(),
                        wind_direction.as_mut_ptr(),
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
    #[test]
    fn atmosphere_new() {
        crate::ceo!(ATMOSPHERE);
    }
}
