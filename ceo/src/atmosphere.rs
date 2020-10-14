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
    built: bool,
}
/// ## `Atmosphere` builder
impl CEO<element::ATMOSPHERE> {
    /// Create a new `Atmosphere` builder
    pub fn new() -> CEO<element::ATMOSPHERE> {
        CEO {
            args: element::ATMOSPHERE {
                r0_at_zenith: 0.16,
                oscale: 25.5,
                zenith_angle: 30f64.to_radians(),
                turbulence: element::TurbulenceProfile::default(),
                ray_tracing: None,
            },
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
            built: true,
        };
        match self.args.ray_tracing {
            None => unsafe {
                atm._c_.setup(
                    self.args.r0_at_zenith as f32,
                    self.args.oscale as f32,
                    self.args.turbulence.n_layer as i32,
                    self.args.turbulence.altitude.as_mut_ptr(),
                    self.args.turbulence.xi0.as_mut_ptr(),
                    self.args.turbulence.wind_speed.as_mut_ptr(),
                    self.args.turbulence.wind_direction.as_mut_ptr(),
                );
            },
            Some(rtc) => match rtc.filepath {
                Some(file) => unsafe {
                    atm._c_.setup2(
                        self.args.r0_at_zenith as f32,
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
                },
                None => unsafe {
                    atm._c_.setup1(
                        self.args.r0_at_zenith as f32,
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
            built: true,
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
        self
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
        self.built = false;
        self
    }
    pub fn gmt_build(&mut self, r_not: f32, l_not: f32) -> &mut Self {
        unsafe {
            self._c_.gmt_setup4(r_not, l_not, 2020);
        }
        self
    }
    /*
    pub fn load_from_json(&mut self, json_file: &str) -> Result<&mut Self, Box<dyn Error>> {
        let mut filename = json_file.to_string();
        filename.push_str(".json");
        let path = Path::new(&filename);
        let file = File::open(path).expect(&format!("{}", path.display()));
        let reader = BufReader::new(file);
        let gmt_atm_args: GmtAtmosphere = serde_json::from_reader(reader)?;
        let ps_path = CString::new(gmt_atm_args.filename.clone()).unwrap();
        unsafe {
            self._c_.gmt_setup6(
                gmt_atm_args.r0,
                gmt_atm_args.l_not,
                gmt_atm_args.length,
                gmt_atm_args.nxy_pupil,
                gmt_atm_args.fov,
                gmt_atm_args.duration,
                ps_path.into_raw(),
                gmt_atm_args.n_duration,
                gmt_atm_args.seed,
            );
        }
        self.filename = String::from(gmt_atm_args.filename);
        self.built = false;
        let ps_path = CString::new(format!("{}", self.filename)).unwrap();
        println!("{:?}", ps_path);
        unsafe {
            self._c_
                .duration_loading(ps_path.into_raw(), self.k_duration);
        }
        Ok(self)
    }
    */
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
    fn drop(&mut self){
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Atmosphere {
    fn time_propagate(&mut self, secs: f64, src: &mut Source) -> &mut Self {
        unsafe {
            let n_xy = src.pupil_sampling;
            let d_xy = (src.pupil_size / (n_xy - 1) as f64) as f32;
            if self.built {
                self._c_
                    .get_phase_screen4(&mut src._c_, d_xy, n_xy, d_xy, n_xy, secs as f32);
            } else {
                /*
                let k_duration = (secs / self._c_.layers_duration as f64) as i32;
                if k_duration > self.k_duration {
                    let ps_path = CString::new(format!("{}", self.filename)).unwrap();
                    //println!("{:?}",ps_path);
                    self._c_.duration_loading(ps_path.into_raw(), k_duration);
                    self.k_duration = k_duration;
                }

                self._c_.ray_tracing(
                    &mut src._c_,
                    d_xy,
                    n_xy,
                    d_xy,
                    n_xy,
                    secs as f32,
                    self.k_duration,
                );
                 */
                self._c_
                    .rayTracing1(&mut src._c_, d_xy, n_xy, d_xy, n_xy, secs as f32);
            }
        }
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

