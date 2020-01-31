use super::{atmosphere, gmt_m1, gmt_m2, modes, shackHartmann, source, vector, zernikeS};

use std::{f32, mem};
use std::ffi::CString;

pub struct OpticsPathSH48 {
    m1_modes: modes,
    m2_modes: zernikeS,
    m1: gmt_m1,
    m2: gmt_m2,
    gs: source,
    wfs: shackHartmann,
    atm: atmosphere,
    pub wfe_rms: f32,
    pupil_size: f64,
    pupil_sampling: i32,
    pub time_stamp: f32,
}
impl OpticsPathSH48 {
    pub fn new() -> OpticsPathSH48 {
        OpticsPathSH48 {
            m1_modes: unsafe { mem::zeroed() },
            m2_modes: unsafe { mem::zeroed() },
            m1: unsafe { mem::zeroed() },
            m2: unsafe { mem::zeroed() },
            gs: unsafe { mem::zeroed() },
            wfs: unsafe { mem::zeroed() },
            atm: unsafe { mem::zeroed() },
            wfe_rms: 0.0,
            pupil_size: 0.0,
            pupil_sampling: 0,
            time_stamp: 0.0,
        }
    }
    pub fn cfg(&mut self) {
        unsafe {
            self.pupil_size = 25.5;
            self.pupil_sampling = 48 * 16 + 1;

            let m1_n_mode = 27;
            let mode_type = CString::new("bending modes").unwrap();
            self.m1_modes
                .setup(mode_type.into_raw() as *mut i8, 7, m1_n_mode);
            self.m1.setup1(&mut self.m1_modes);

            let mut a = vec![0.0];
            self.m2_modes.setup1(0, a.as_mut_ptr(), 7);
            self.m2.setup2(&mut self.m2_modes);

            let origin = vector {
                x: 0.0,
                y: 0.0,
                z: 25.0,
            };
            let band = String::from("V");
            let mut magnitude: Vec<f32> = vec![0.0];
            let mut zenith: Vec<f32> = vec![0.0];
            let mut azimuth: Vec<f32> = vec![0.0];
            self.gs.setup7(
                band.as_ptr() as *const i8,
                magnitude.as_mut_ptr(),
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                f32::INFINITY,
                zenith.len() as i32,
                self.pupil_size,
                self.pupil_sampling,
                origin,
            );

            let d = self.pupil_size / 48.0;
            self.wfs.setup(48, 16, d as f32, 2, 24, 3, 3);

            self.goxp();
            self.gs.fwhm = 3.16;
            self.wfs.calibrate(&mut self.gs, 0.0);
            self.gs.fwhm = 0.0;

            self.atm.gmt_setup6(
                15e-2,
                25.0,
                26.0,
                351,
                20.0 * f32::consts::PI / 180.0 / 60.0,
                15.0,
                String::from("/home/rconan/DATA/gmtAtmosphereL030_2.bin").as_mut_ptr() as *mut i8,
                4,
                2019,
            );
        }
    }
    pub fn goxp(&mut self) {
        unsafe {
            let d = (self.pupil_size / ((self.pupil_sampling - 1) as f64)) as f32;
            self.gs.wavefront.reset();
            self.gs.reset_rays();
            self.m2.blocking(&mut self.gs.rays);
            self.m1.trace(&mut self.gs.rays);
            self.gs.rays.gmt_truss_onaxis();
            self.m2.trace(&mut self.gs.rays);
            self.gs.rays.to_sphere1(-5.830, 2.197173);
            self.gs.opd2phase();
            self.atm.rayTracing1(
                &mut self.gs,
                d,
                self.pupil_sampling,
                d,
                self.pupil_sampling,
                self.time_stamp,
            );
            self.gs.wavefront.rms(&mut self.wfe_rms);
        }
    }
    pub fn gofp(&mut self) {
        self.goxp();
        unsafe { self.wfs.propagate(&mut self.gs) }
    }
}
impl Drop for OpticsPathSH48 {
    fn drop(&mut self) {
        unsafe {
            self.m1_modes.cleanup();
            self.m1.cleanup();
            self.m2_modes.cleanup();
            self.m2.cleanup();
            self.gs.cleanup();
            self.wfs.cleanup();
        }
    }
}
