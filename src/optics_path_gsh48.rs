use super::{
    geometricShackHartmann, gmt_m1, gmt_m2, modes, source, vector, zernikeS,dev2host
};

use std::{f32, mem};

pub struct OpticsPathGSH48 {
    m1_modes: modes,
    m2_modes: zernikeS,
    pub m1: gmt_m1,
    pub m2: gmt_m2,
    gs: source,
    wfs: geometricShackHartmann,
    pub wfe_rms: f32,
}
impl OpticsPathGSH48 {
    pub fn new() -> OpticsPathGSH48 {
        OpticsPathGSH48 {
            m1_modes: unsafe { mem::zeroed() },
            m2_modes: unsafe { mem::zeroed() },
            m1: unsafe { mem::zeroed() },
            m2: unsafe { mem::zeroed() },
            gs: unsafe { mem::zeroed() },
            wfs: unsafe { mem::zeroed() },
            wfe_rms: 0.0,
        }
    }
    pub fn cfg(&mut self, n_side_lenslet: i32, n_px_lenslet: i32) {
        unsafe {
            let pupil_size = 25.5;
            let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;

            let m1_n_mode = 27;
            let mut mode_type = String::from("bending modes");
            self.m1_modes
                .setup(mode_type.as_mut_ptr() as *mut i8, 7, m1_n_mode);
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
                pupil_size,
                pupil_sampling,
                origin,
            );

            let d = pupil_size / 48.0;
            self.wfs.setup(48, d as f32, 3);

            self.go_xp();
            self.gs.fwhm = 3.16;
            self.wfs.calibrate(&mut self.gs, 0.0);
            self.gs.fwhm = 0.0;
        }
    }
    pub fn go_xp(&mut self) {
        unsafe {
            self.gs.wavefront.reset();
            self.gs.reset_rays();
            self.m2.blocking(&mut self.gs.rays);
            self.m1.trace(&mut self.gs.rays);
            self.gs.rays.gmt_truss_onaxis();
            self.m2.trace(&mut self.gs.rays);
            self.gs.rays.to_sphere1(-5.830, 2.197173);
            self.gs.opd2phase();
        };
    }
    pub fn go_fp(&mut self) {
        self.go_xp();
        unsafe {
            self.wfs.propagate(&mut self.gs);
            self.wfs.process();
            self.wfs.camera.reset();
        }
    }
    pub fn m1_status(&mut self) {
        unsafe { self.m1.motion_CS.info_details() };
    }
    pub fn m1_update(&mut self, o: [f64; 3], a: [f64; 3], seg_id: i32) {
        let t_xyz = vector {
            x: o[0],
            y: o[1],
            z: o[2],
        };
        let r_xyz = vector {
            x: a[0],
            y: a[1],
            z: a[2],
        };
        unsafe { self.m1.motion_CS.update1(t_xyz, r_xyz, seg_id); }
    }
    pub fn centroids(&mut self, c: &mut Vec<f32>, n_c: i32) {
        unsafe {
            dev2host(c.as_mut_ptr(),self.wfs.data_proc.d__c,n_c);
        }
    }
}
impl Drop for OpticsPathGSH48 {
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
