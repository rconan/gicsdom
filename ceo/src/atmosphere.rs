use serde::Deserialize;
use std::{f32, mem};

use super::ceo_bindings::atmosphere;
use super::Propagation;
use super::Source;

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
    ) -> &mut Self {
        unsafe {
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
impl Propagation for Atmosphere {
    
    fn time_propagate(&mut self, secs: f64, src: &mut Source) -> &mut Self {
        unsafe {
            //src._c_.wavefront.reset();
            let n_xy = src.pupil_sampling;
            let d_xy = (src.pupil_size / (n_xy - 1) as f64) as f32;
            if self.built {
                self._c_
                    .get_phase_screen4(&mut src._c_, d_xy, n_xy, d_xy, n_xy, 0f32);
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
            }
            //src._c_.opd2phase();
        }
        self
    }
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        self.time_propagate(self.secs, src)
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ceo_built_atmosphere() {
        let mut gmt = Gmt::new(0, None);
        gmt.build();
        let mut wfs = ShackHartmann::new(1, 48, 16, 25.5 / 48.0);
        wfs.build(8, Some(24), None);
        let mut src = wfs.new_guide_stars();
        src.build("V", vec![0.0f32], vec![0.0f32], vec![0.0f32]);
        let mut atm = Atmosphere::new();
        atm.build(0.15, 25.);
        let n = 10;
        let mut wfe_rms = (0..n)
            .into_iter()
            .map(|i| {
                atm.secs = i as f64;
                src.through(&mut gmt).xpupil().through(&mut atm);
                src.wfe_rms_10e(-9)[0]
            })
            .collect::<Vec<f32>>();
        wfe_rms.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("WFE RMS: {:?}nm", wfe_rms);
    }
    #[test]
    fn ceo_load_atmosphere() {
        let mut gmt = Gmt::new(0, None);
        gmt.build();
        let mut wfs = ShackHartmann::new(1, 48, 16, 25.5 / 48.0);
        wfs.build(8, Some(24), None);
        let mut src = wfs.new_guide_stars();
        src.build("V", vec![0.0f32], vec![0.0f32], vec![0.0f32]);
        let mut atm = Atmosphere::new();
        atm.load_from_json("/home/ubuntu/DATA/gmtAtmosphereL025_1579821046.json")
            .unwrap();
        let n = 10;
        let mut wfe_rms = (0..n)
            .into_iter()
            .map(|i| {
                atm.secs = i as f64;
                src.through(&mut gmt).xpupil().through(&mut atm);
                src.wfe_rms_10e(-9)[0]
            })
            .collect::<Vec<f32>>();
        wfe_rms.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("WFE RMS: {:?}nm", wfe_rms);
    }
}
*/
