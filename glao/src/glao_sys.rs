use crate::system::{GlaoField, System};
use ceo::{
    pssn::AtmosphereTelescopeError as ATE, pssn::TelescopeError as TE, Atmosphere, Conversion, Cu,
    Gmt, Mask, PSSn, Source,
};
use cfd::DomeSeeing;
use serde::ser::{Serialize, SerializeStruct, Serializer};
use serde::{Deserialize as Dserde, Serialize as Sserde};
use serde_pickle as pickle;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;
use std::{f32, f64};

const PUPIL_SIZE: f64 = 25.5;

#[derive(Sserde, Dserde, Default, Debug)]
pub struct Fwhms {
    atm_fwhm_x: f64,
    atm_fwhm: Vec<f64>,
    glao_fwhm: Vec<f64>,
}

pub struct ScienceField {
    band: String,
    n_src: usize,
    n_px: usize,
    zenith: Vec<f32>,
    azimuth: Vec<f32>,
    pub src: Source,
    atm_pssn: ceo::PSSn<ATE>,
    pub pssn: ceo::PSSn<ATE>,
    pub pssn_nsample_tol: Option<(usize, f32)>,
    fwhms: Fwhms,
}
impl ScienceField {
    pub fn on_axis(band: &str, n_px: usize, r0_oscale: Option<(f32, f32)>) -> Self {
        ScienceField {
            band: band.to_owned(),
            n_src: 1,
            n_px,
            zenith: vec![0f32],
            azimuth: vec![0f32],
            src: Source::new(1, PUPIL_SIZE, n_px as i32),
            atm_pssn: match r0_oscale {
                Some(r0_oscale) => PSSn::from_r0_and_outerscale(r0_oscale.0, r0_oscale.1),
                None => PSSn::new(),
            },
            pssn: match r0_oscale {
                Some(r0_oscale) => PSSn::from_r0_and_outerscale(r0_oscale.0, r0_oscale.1),
                None => PSSn::new(),
            },
            pssn_nsample_tol: None,
            fwhms: Fwhms::default(),
        }
    }
    pub fn delaunay_21(band: &str, n_px: usize, r0_oscale: Option<(f32, f32)>) -> Self {
        let glao_field_reader = File::open("glao_field.pkl").expect("File not found!");
        let glao_field: GlaoField =
            pickle::from_reader(glao_field_reader).expect("File loading failed!");
        let n_src = glao_field.zenith_arcmin.len();
        assert_eq!(n_src, 21);
        let src_zen = glao_field
            .zenith_arcmin
            .iter()
            .map(|x| x.from_arcmin())
            .collect::<Vec<f32>>();
        let src_azi = glao_field
            .azimuth_degree
            .iter()
            .map(|x| x.to_radians())
            .collect::<Vec<f32>>();
        ScienceField {
            band: band.to_owned(),
            n_src: 21,
            n_px,
            zenith: src_zen,
            azimuth: src_azi,
            src: Source::new(n_src as i32, PUPIL_SIZE, n_px as i32),
            atm_pssn: match r0_oscale {
                Some(r0_oscale) => PSSn::from_r0_and_outerscale(r0_oscale.0, r0_oscale.1),
                None => PSSn::new(),
            },
            pssn: match r0_oscale {
                Some(r0_oscale) => PSSn::from_r0_and_outerscale(r0_oscale.0, r0_oscale.1),
                None => PSSn::new(),
            },
            pssn_nsample_tol: None,
            fwhms: Fwhms::default(),
        }
    }
    pub fn build(&mut self) -> &mut Self {
        self.src.build(
            &self.band,
            self.zenith.clone(),
            self.azimuth.clone(),
            vec![0f32; self.n_src],
        );
        let mut gmt = Gmt::new();
        gmt.build(0, Some(0));
        self.src.through(&mut gmt).xpupil();
        self.pssn.build(&mut self.src);
        self.atm_pssn.build(&mut self.src);
        self
    }
    pub fn wrap_up(&mut self) -> &mut Self {
        self.atm_pssn.peek().telescope_error_into_otf();
        self.pssn.peek().telescope_error_into_otf();
        let mut fwhm = ceo::Fwhm::new();
        fwhm.build(&mut self.src);
        fwhm.upper_bracket = 2f64 / self.pssn.r0() as f64;
        let atm_fwhm_x = ceo::Fwhm::atmosphere(
            self.pssn.wavelength as f64,
            self.pssn.r0() as f64,
            self.pssn.oscale as f64,
        )
        .to_arcsec();
        let atm_fwhm = fwhm.from_complex_otf(&self.atm_pssn.telescope_error_otf());
        let glao_fwhm = fwhm.from_complex_otf(&self.pssn.telescope_error_otf());
        self.fwhms.atm_fwhm_x = atm_fwhm_x;
        self.fwhms.atm_fwhm = atm_fwhm.iter().map(|x| x.to_arcsec()).collect::<Vec<_>>();
        self.fwhms.glao_fwhm = glao_fwhm.iter().map(|x| x.to_arcsec()).collect::<Vec<_>>();
        self
    }
    pub fn glao_wrap_up(&mut self) -> (f64,Vec<f64>) {
        let mut fwhm = ceo::Fwhm::new();
        fwhm.build(&mut self.src);
        fwhm.upper_bracket = 2f64 / self.pssn.r0() as f64;
        let atm_fwhm_x = ceo::Fwhm::atmosphere(
            self.pssn.wavelength as f64,
            self.pssn.r0() as f64,
            self.pssn.oscale as f64,
        )
            .to_arcsec();
        let fwhm = fwhm.from_complex_otf(&self.pssn.telescope_error_otf());
        (atm_fwhm_x,fwhm.iter().map(|x| x.to_arcsec()).collect::<Vec<_>>())
    }
    pub fn results(&self) -> (Vec<f32>, Vec<f32>, f64, Vec<f64>, Vec<f64>) {
        (
            self.atm_pssn.estimates.clone(),
            self.pssn.estimates.clone(),
            self.fwhms.atm_fwhm_x,
            self.fwhms.atm_fwhm.clone(),
            self.fwhms.glao_fwhm.clone(),
        )
    }
    pub fn dump(&self, filename: &str) -> Result<&Self, Box<dyn Error>> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::with_capacity(1_000_000, file);
        pickle::to_writer(&mut writer, self, true)?;
        Ok(self)
    }
}

impl Serialize for ScienceField {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("ScienceField", 5)?;
        state.serialize_field("pupil sampling [px]", &self.n_px)?;
        state.serialize_field("atmosphere pssn", &self.atm_pssn)?;
        state.serialize_field("GLAO pssn", &self.pssn)?;
        state.serialize_field("FWHMS [arcsec]", &self.fwhms)?;
        state.serialize_field("pssn (n_sample,tol)", &self.pssn_nsample_tol)?;
        state.end()
    }
}

pub struct GlaoSys<'a, 'b> {
    pub sys: System,
    lenslet_mask: Mask,
    calib: Cu<f32>,
    step: usize,
    pub atm: &'a mut Atmosphere,
    pub science: &'b mut ScienceField,
    d_mean_c: ceo::Cu<f32>,
    x: ceo::Cu<f32>,
    pub s12: Option<((usize, usize), (usize, usize))>,
}

impl<'a, 'b> GlaoSys<'a, 'b> {
    pub fn new(
        pupil_size: f64,
        n_sensor: i32,
        n_lenslet: i32,
        n_px_lenslet: i32,
        science: &'b mut ScienceField,
        atm: &'a mut Atmosphere,
    ) -> Self {
        Self {
            sys: System::new(pupil_size, n_sensor as i32, n_lenslet, n_px_lenslet),
            lenslet_mask: Mask::new(),
            calib: Cu::new(),
            step: 0,
            atm,
            science,
            d_mean_c: Cu::new(),
            x: Cu::new(),
            s12: None,
        }
    }
    pub fn default(atm: &'a mut Atmosphere, science: &'b mut ScienceField) -> Self {
        Self {
            sys: System::new(25.5, 4, 48, 16),
            lenslet_mask: Mask::new(),
            calib: Cu::new(),
            step: 0,
            atm,
            science,
            d_mean_c: Cu::new(),
            x: Cu::new(),
            s12: None,
        }
    }
    pub fn build(
        &mut self,
        gs_off_axis_dist: f32,
        m2_n_kl: usize,
        wfs_intensity_threshold: f64,
    ) -> &mut Self {
        let n_gs = self.sys.n_wfs;
        let gs_zen = (0..n_gs).map(|_| gs_off_axis_dist).collect::<Vec<f32>>();
        let gs_azi = (0..n_gs)
            .map(|x| (x as f32) * 2f32 * f32::consts::PI / n_gs as f32)
            .collect::<Vec<f32>>();
        self.sys
            .gmt_build("S12", 2, m2_n_kl)
            .wfs_build("Vs", gs_zen, gs_azi, vec![0f32; n_gs])
            .wfs_calibrate(wfs_intensity_threshold)
            .through();
        log::info!("GLAO centroids #: {}", self.sys.wfs.n_centroids);

        let n_lenslet = self.sys.n_lenslet as usize;
        let (_, _, wfs) = self.sys.devices();
        let m = n_lenslet * n_lenslet;
        let mask = wfs
            .lenset_mask()
            .from_dev()
            .chunks(m)
            .fold(vec![0u8; m], |a, x| {
                a.iter()
                    .zip(x.iter())
                    .map(|y: (&u8, &f32)| if *y.1 > 0f32 { 1u8 } else { *y.0 })
                    .collect::<Vec<u8>>()
            });
        let nnz = mask.iter().cloned().map(|x| x as usize).sum::<usize>();
        log::info!("Centroid mask: [{}], nnz: {}", mask.len(), nnz);

        self.lenslet_mask.build(m);
        let mut f: ceo::Cu<f32> = ceo::Cu::vector(m);
        let mut q = mask.iter().cloned().map(|x| x as f32).collect::<Vec<f32>>();
        f.to_dev(&mut q);
        self.lenslet_mask.filter(&mut f);

        self
    }
    pub fn calibration(&mut self) -> &mut Self {
        let now = Instant::now();
        self.calib = {
            let mut on_axis_sys = System::new(
                self.sys.pupil_size,
                1,
                self.sys.n_lenslet as i32,
                self.sys.n_px_lenslet as i32,
            );
            on_axis_sys
                .gmt_build("bending modes", 27, self.sys.m2_n_mode)
                .wfs_build("Vs", vec![0f32], vec![0f32], vec![0f32])
                .wfs_calibrate(0f64)
                .through();
            on_axis_sys.m2_mode_calibrate_data(&mut self.lenslet_mask)
        };
        self.calib.qr();
        log::info!("Calibration & Inversion in {}s", now.elapsed().as_secs());
        log::info!(
            "Calibration matrix size [{}x{}]",
            self.calib.n_row(),
            self.calib.n_col()
        );
        self.d_mean_c = Cu::vector(self.calib.n_row());
        self.d_mean_c.malloc();
        self.x = Cu::vector(self.calib.n_col());
        self.x.malloc();
        self
    }
    pub fn peek(&mut self) -> Vec<f32> {
        self.science.pssn.peek().estimates.clone()
    }
    pub fn set_m1_polishing_error(&mut self, v: f64) {
        if let Some(((s1, e1), (s2, e2))) = self.s12 {
            let mut a = vec![vec![0f64; self.sys.gmt.m1_n_mode]; 7];
            ((s1 - 1)..e1).for_each(|i| {
                a[i][0] = v;
            });
            ((s2 - 1)..e2).for_each(|i| {
                a[i][1] = v;
            });
        }
    }
    pub fn get_science(&mut self, disturbances: &mut DomeSeeing) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
        let (_gs, gmt, _wfs) = self.sys.devices();
        let wfe_rms = self
            .science
            .src
            .through(gmt)
            .xpupil()
            .through(self.atm)
            .through(disturbances)
            .through(&mut self.science.pssn)
            .wfe_rms_10e(-9);
        (
            wfe_rms,
            self.science.src.segment_wfe_rms_10e(-9),
            self.science.src.segment_piston_10e(-9),
        )
    }
    pub fn closed_loop(
        &mut self,
        disturbances: &mut DomeSeeing,
        kl_coefs: &mut Vec<Vec<f64>>,
        gain: f64,
    ) {
        let (gs, gmt, wfs) = self.sys.devices();

        wfs.reset();
        gs.through(gmt)
            .xpupil()
            .through(self.atm)
            .through(disturbances)
            .through(wfs);
        wfs.process()
            .fold_into(&mut self.d_mean_c, &mut self.lenslet_mask);

        self.calib.qr_solve_as_ptr(&mut self.x, &mut self.d_mean_c);
        let h_x = self.x.from_dev();
        let mut k = 0;
        for s in kl_coefs.iter_mut() {
            for a in s.iter_mut().skip(1) {
                *a -= gain * h_x[k] as f64;
                k += 1;
            }
        }

        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);
    }
}

pub fn m1_polishing_wavefront_error(glao: &mut GlaoSys) {
    glao.set_m1_polishing_error(1f64);
    let (_, gmt, _) = glao.sys.devices();
    let mut src = Source::new(1, 25.5, 1024);
    src.build("Vs", vec![0f32], vec![0f32], vec![0f32]);
    let wfe_rms = src.through(gmt).xpupil().wfe_rms_10e(-9);
    log::info!("M1 polishing WFE RMS [nm]: {:.3}", wfe_rms[0]);
    let mut file = File::create("m1_polishing_wavefront_error.pkl").unwrap();
    pickle::to_writer(&mut file, src.phase(), true).unwrap();
}

impl<'a, 'b> Iterator for GlaoSys<'a, 'b> {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        let n_kl = self.sys.m2_n_mode;

        self.set_m1_polishing_error(0f64);

        let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        self.sys.gmt.set_m2_modes(&mut b);

        self.science
            .src
            .through(&mut self.sys.gmt)
            .xpupil()
            .through(self.atm)
            .through(&mut self.science.atm_pssn);

        self.set_m1_polishing_error(1f64);
        let (gs, gmt, wfs) = self.sys.devices();

        wfs.reset();
        gs.through(gmt).xpupil().through(self.atm).through(wfs);
        wfs.process()
            .fold_into(&mut self.d_mean_c, &mut self.lenslet_mask);

        self.calib.qr_solve_as_ptr(&mut self.x, &mut self.d_mean_c);
        let h_x = self.x.from_dev();
        let mut k = 0;
        for s in kl_coefs.iter_mut() {
            for a in s.iter_mut().skip(1) {
                *a -= h_x[k] as f64;
                k += 1;
            }
        }

        let mut b = kl_coefs.into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);

        self.science
            .src
            .through(gmt)
            .xpupil()
            .through(self.atm)
            .through(&mut self.science.pssn);

        /*let wfe_rms = src.wfe_rms_10e(-9)[0];
        if k_sample % 100 == 0 {
            println!(
                "#{:6}: WFE RMS: {:5.0}/{:5.0}nm ; PSSn: {}",
                k_sample,
                wfe_rms_0,
                wfe_rms,
                pssn.peek()
            );
        }*/

        self.atm.reset();

        self.step += 1;
        Some(self.step)
    }
}
