use super::bindings::*;

use std::ffi::CString;
use std::{f32, mem};

pub struct Gmt {
    _c_m1_modes: modes,
    _c_m2_modes: zernikeS,
    _c_m1: gmt_m1,
    _c_m2: gmt_m2,
    m1_n_mode: u64,
    m2_n_mode: u64,
    a: Vec<f64>,
}
impl Gmt {
    pub fn test(m1_n_mode: u64, m2_n_mode: Option<u64>) -> Gmt {
        let mode_type = CString::new("bending modes").unwrap();
        let mut this: Gmt = Gmt {
            _c_m1_modes: unsafe { mem::zeroed() },
            _c_m2_modes: unsafe { mem::zeroed() },
            _c_m1: unsafe { mem::zeroed() },
            _c_m2: unsafe { mem::zeroed() },
            m1_n_mode,
            m2_n_mode: match m2_n_mode {
                Some(m2_n_mode) => m2_n_mode,
                None => 0,
            },
            a: vec![0.0],
        };
        if this.m2_n_mode > 0 {
            this.a = vec![0.0; this.m2_n_mode as usize];
        }
        unsafe {
            this._c_m1_modes
                .setup(mode_type.into_raw(), 7, this.m1_n_mode as i32);
            this._c_m1.setup1(&mut this._c_m1_modes);
            this._c_m2_modes
                .setup1(this.m2_n_mode as i32, this.a.as_mut_ptr(), 7);
            this._c_m2.setup2(&mut this._c_m2_modes);
        }
        this
    }
    pub fn new(m1_n_mode: u64, m2_n_mode: Option<u64>) -> Gmt {
        Gmt {
            _c_m1_modes: unsafe { mem::zeroed() },
            _c_m2_modes: unsafe { mem::zeroed() },
            _c_m1: unsafe { mem::zeroed() },
            _c_m2: unsafe { mem::zeroed() },
            m1_n_mode,
            m2_n_mode: match m2_n_mode {
                Some(m2_n_mode) => m2_n_mode,
                None => 0,
            },
            a: vec![0.0],
        }
    }
    pub fn build(&mut self) -> &mut Gmt {
        let mode_type = CString::new("bending modes").unwrap();
        if self.m2_n_mode > 0 {
            self.a = vec![0.0; self.m2_n_mode as usize];
        }
        unsafe {
            self._c_m1_modes
                .setup(mode_type.into_raw(), 7, self.m1_n_mode as i32);
            self._c_m1.setup1(&mut self._c_m1_modes);
            self._c_m2_modes
                .setup1(self.m2_n_mode as i32, self.a.as_mut_ptr(), 7);
            self._c_m2.setup2(&mut self._c_m2_modes);
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        let mut a: Vec<f64> = vec![0.0; 7 * self.m1_n_mode as usize];
        unsafe {
            self._c_m1.reset();
            self._c_m2.reset();
            self._c_m1_modes.update(a.as_mut_ptr());
        }
        self
    }
    pub fn set_m1_segment_state(&mut self, sid: i32, t_xyz: &Vec<f64>, r_xyz: &Vec<f64>) {
        let t_xyz = vector {
            x: t_xyz[0],
            y: t_xyz[1],
            z: t_xyz[2],
        };
        let r_xyz = vector {
            x: r_xyz[0],
            y: r_xyz[1],
            z: r_xyz[2],
        };
        unsafe {
            self._c_m1.update(t_xyz, r_xyz, sid);
        }
    }
    pub fn set_m2_segment_state(&mut self, sid: i32, t_xyz: &Vec<f64>, r_xyz: &Vec<f64>) {
        let t_xyz = vector {
            x: t_xyz[0],
            y: t_xyz[1],
            z: t_xyz[2],
        };
        let r_xyz = vector {
            x: r_xyz[0],
            y: r_xyz[1],
            z: r_xyz[2],
        };
        unsafe {
            self._c_m2.update(t_xyz, r_xyz, sid);
        }
    }
    pub fn set_m1_modes(&mut self, a: &mut Vec<f64>) {
        unsafe {
            self._c_m1_modes.update(a.as_mut_ptr());
        }
    }
}
impl Drop for Gmt {
    fn drop(&mut self) {
        unsafe {
            self._c_m1_modes.cleanup();
            self._c_m1.cleanup();
            self._c_m2_modes.cleanup();
            self._c_m2.cleanup();
        }
    }
}
pub struct Source {
    _c_: source,
    size: i32,
    pupil_size: f64,
    pupil_sampling: i32,
    pub _wfe_rms: Vec<f32>,
    is_ray_trace: bool,
    is_opd_to_phase: bool,
}
impl Source {
    pub fn test(
        band: &str,
        zenith: &mut Vec<f32>,
        azimuth: &mut Vec<f32>,
        magnitude: &mut Vec<f32>,
        pupil_size: f64,
        pupil_sampling: i32,
    ) -> Source {
        assert_eq!(zenith.len(), azimuth.len());
        assert_eq!(zenith.len(), magnitude.len());
        let band = CString::new(band).unwrap();
        let size = zenith.len() as i32;
        let mut this = Source {
            _c_: unsafe { mem::zeroed() },
            size: size,
            pupil_size,
            pupil_sampling,
            _wfe_rms: Vec::new(),
            is_ray_trace: false,
            is_opd_to_phase: false,
        };
        let mut_ref_to_this = &mut this;
        unsafe {
            let origin = vector {
                x: 0.0,
                y: 0.0,
                z: 25.0,
            };
            mut_ref_to_this._c_.setup7(
                band.into_raw(),
                magnitude.as_mut_ptr(),
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                f32::INFINITY,
                size,
                pupil_size,
                pupil_sampling,
                origin,
            );
        }
        this
    }
    pub fn new(size: i32, pupil_size: f64, pupil_sampling: i32) -> Source {
        Source {
            _c_: unsafe { mem::zeroed() },
            size,
            pupil_size,
            pupil_sampling,
            _wfe_rms: vec![0.0; size as usize],
            is_ray_trace: false,
            is_opd_to_phase: false,
        }
    }
    pub fn from(args: (i32, f64, i32)) -> Source {
        Source::new(args.0, args.1, args.2)
    }
    pub fn build(
        &mut self,
        band: &str,
        mut zenith: Vec<f32>,
        mut azimuth: Vec<f32>,
        mut magnitude: Vec<f32>,
    ) -> &mut Self {
        assert_eq!(zenith.len(), azimuth.len());
        assert_eq!(zenith.len(), magnitude.len());
        let band = CString::new(band).unwrap();
        unsafe {
            let origin = vector {
                x: 0.0,
                y: 0.0,
                z: 25.0,
            };
            self._c_.setup7(
                band.into_raw(),
                magnitude.as_mut_ptr(),
                zenith.as_mut_ptr(),
                azimuth.as_mut_ptr(),
                f32::INFINITY,
                self.size,
                self.pupil_size,
                self.pupil_sampling,
                origin,
            );
        }
        self
    }
    pub fn wfe_rms(&mut self) {
        if self.is_ray_trace {
            if !self.is_opd_to_phase {
                unsafe {
                    self._c_.wavefront.reset();
                    self._c_.opd2phase();
                }
            }
            unsafe {
                self._c_.wavefront.rms(self._wfe_rms.as_mut_ptr());
            }
        } else {
            println!("Ray tracing through the GMT must be performed first!")
        }
    }
    pub fn wfe_rms_10e(&mut self, exp: i32) -> Vec<f32> {
        self.wfe_rms();
        self._wfe_rms
            .iter()
            .map(|x| x * 10_f32.powi(-exp))
            .collect()
    }
    pub fn reset(&mut self) {
        unsafe {
            self._c_.wavefront.reset();
            self._c_.reset_rays();
        }
        self.is_ray_trace = false;
        self.is_opd_to_phase = false;
    }
}
impl Drop for Source {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
pub trait Propagation {
    fn propagate(&mut self, src: &mut Source) -> &mut Self;
}
impl Propagation for Gmt {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            src._c_.reset_rays();
            let rays: &mut bundle = &mut src._c_.rays;
            self._c_m2.blocking(rays);
            self._c_m1.trace(rays);
            //gs.rays.gmt_truss_onaxis();
            self._c_m2.trace(rays);
            rays.to_sphere1(-5.830, 2.197173);
        }
        src.is_ray_trace = true;
        src.is_opd_to_phase = false;
        self
    }
}
impl Source {
    pub fn through<T: Propagation>(&mut self, system: &mut T) -> &mut Self {
        system.propagate(self);
        self
    }
}
pub struct PSSn {
    _c_: pssn,
    r_not: f64,
    l_not: f64,
    zenith: f64,
}
impl PSSn {
    pub fn new(r_not: f64, l_not: f64, zenith: f64) -> PSSn {
        PSSn {
            _c_: unsafe{ mem::zeroed() },
            r_not,
            l_not,
            zenith,
        }
    }
    pub fn build(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.setup(&mut src._c_, self.r_not as f32, self.l_not as  f32);
        }
        self
    }
    pub fn reset(&mut self, src: &mut Source) -> f32 {
        let pssn_val: f32;
        unsafe {
            self._c_.otf(&mut src._c_);
            pssn_val = self._c_.eval();
            self._c_.N_O = 0;
        }
        pssn_val
    }
    pub fn peek(&mut self, src: &mut Source) -> f32 {
        let pssn_val: f32;
        unsafe {
            self._c_.otf(&mut src._c_);
            pssn_val = self._c_.eval();
        }
        pssn_val
    }
    pub fn sum(&mut self, src: &mut Source) {
        unsafe {
            self._c_.otf(&mut src._c_);
        }
    }
}
pub struct Geometric_ShackHartmann {
    _c_: geometricShackHartmann,
    n_side_lenslet: i32,
    n_px_lenslet: i32,
    d: f64,
    n_sensor: i32,
    pub n_centroids: i32,
    pub centroids: Vec<f32>,
}
pub struct Diffractive_ShackHartmann {
    _c_: shackHartmann,
    n_side_lenslet: i32,
    n_px_lenslet: i32,
    d: f64,
    n_sensor: i32,
}
impl Geometric_ShackHartmann {
    pub fn new(
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        d: f64,
    ) -> Geometric_ShackHartmann {
        Geometric_ShackHartmann {
            _c_: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
            n_centroids: 0,
            centroids: Vec::new(),
        }
    }
    pub fn build(&mut self) -> &mut Self {
        unsafe {
            self._c_
                .setup(self.n_side_lenslet, self.d as f32, self.n_sensor);
        }
        self.n_centroids = self.n_side_lenslet * self.n_side_lenslet * 2 * self.n_sensor;
        self.centroids = vec![0.0; self.n_centroids as usize]; //Vec::with_capacity(self.n_centroids as usize);
        self
    }
    pub fn guide_star_args(&self) -> (i32, f64, i32) {
        (
            self.n_sensor,
            self.d * self.n_side_lenslet as f64,
            self.n_px_lenslet * self.n_side_lenslet + 1,
        )
    }
    pub fn new_guide_stars(&self) -> Source {
        Source::new(
            self.n_sensor,
            self.d * self.n_side_lenslet as f64,
            self.n_px_lenslet * self.n_side_lenslet + 1,
        )
    }
    pub fn calibrate(&mut self, src: &mut Source, threshold: f64) -> Result<&mut Self, String> {
        if src.is_ray_trace {
            if !src.is_opd_to_phase {
                unsafe {
                    src._c_.wavefront.reset();
                    src._c_.opd2phase();
                }
            }
            unsafe {
                self._c_.calibrate(&mut src._c_, threshold as f32);
            }
        } else {
            return Err("Ray tracing through the GMT must be performed first!".to_string());
        }
        Ok(self)
    }
    pub fn process(&mut self) -> &mut Self {
        self.centroids = vec![0.0; self.n_centroids as usize];
        //let mut c = vec![0.0;self.n_centroids as usize];
        unsafe {
            self._c_.process();
            dev2host(
                self.centroids.as_mut_ptr(),
                self._c_.data_proc.d__c,
                self.n_centroids as i32,
            );
            self._c_.reset();
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.reset();
        }
        self
    }
}
impl Drop for Geometric_ShackHartmann {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Geometric_ShackHartmann {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        if src.is_ray_trace {
            if !src.is_opd_to_phase {
                unsafe {
                    src._c_.wavefront.reset();
                    src._c_.opd2phase();
                }
            }
            unsafe {
                self._c_.propagate(&mut src._c_);
            }
        } else {
            println!("Ray tracing through the GMT must be performed first!")
        }
        self
    }
}
impl Diffractive_ShackHartmann {
    pub fn new(
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        d: f64,
    ) -> Diffractive_ShackHartmann {
        Diffractive_ShackHartmann {
            _c_: unsafe { mem::zeroed() },
            n_side_lenslet,
            n_px_lenslet,
            d,
            n_sensor,
        }
    }
    pub fn build(
        &mut self,
        n_px_framelet: i32,
        n_px_imagelet: Option<i32>,
        osf: Option<i32>,
    ) -> &mut Self {
        let n_px = match n_px_imagelet {
            Some(n_px_imagelet) => n_px_imagelet,
            None => n_px_framelet,
        };
        let b = n_px / n_px_framelet;
        let o = match osf {
            Some(osf) => osf,
            None => 2,
        };
        unsafe {
            self._c_.setup(
                self.n_side_lenslet,
                self.n_px_lenslet,
                self.d as f32,
                o,
                n_px,
                b,
                self.n_sensor,
            );
        }
        self
    }
    pub fn guide_star_args(&self) -> (i32, f64, i32) {
        (
            self.n_sensor,
            self.d * self.n_side_lenslet as f64,
            self.n_px_lenslet * self.n_side_lenslet + 1,
        )
    }
    pub fn calibrate(&mut self, src: &mut Source, threshold: f64) -> &mut Self {
        if src.is_ray_trace {
            if !src.is_opd_to_phase {
                unsafe {
                    src._c_.wavefront.reset();
                    src._c_.opd2phase();
                }
            }
            unsafe {
                self._c_.calibrate(&mut src._c_, threshold as f32);
                self._c_.camera.reset();
            }
        } else {
            println!("Ray tracing through the GMT must be performed first!")
        }
        self
    }
    pub fn process(&mut self) -> &mut Self {
        unsafe {
            self._c_.process();
        }
        self
    }
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.camera.reset();
        }
        self
    }
}
impl Drop for Diffractive_ShackHartmann {
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Diffractive_ShackHartmann {
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        if src.is_ray_trace {
            if !src.is_opd_to_phase {
                unsafe {
                    src._c_.wavefront.reset();
                    src._c_.opd2phase();
                }
            }
            unsafe {
                self._c_.propagate(&mut src._c_);
            }
        } else {
            println!("Ray tracing through the GMT must be performed first!")
        }
        self
    }
}
