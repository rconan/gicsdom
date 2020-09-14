use crate::system::System;
use ceo::{pssn::AtmosphereTelescopeError as ATE, Atmosphere, Conversion, Cu, Mask, PSSn, Source};
use std::f32;
use std::time::Instant;

pub struct GlaoSys<'a> {
    sys: System,
    lenslet_mask: Mask,
    calib: Cu<f32>,
    src: Source,
    atm_pssn: ceo::PSSn<ATE>,
    pssn: ceo::PSSn<ATE>,
    step: usize,
    atm: &'a mut Atmosphere,
    d_mean_c: ceo::Cu<f32>,
    x: ceo::Cu<f32>,
}

impl<'a> GlaoSys<'a> {
    pub fn new(
        pupil_size: f64,
        n_sensor: i32,
        n_lenslet: i32,
        n_px_lenslet: i32,
        n_src: usize,
        atm: &'a mut Atmosphere,
    ) -> Self {
        Self {
            sys: System::new(pupil_size, n_sensor as i32, n_lenslet, n_px_lenslet),
            lenslet_mask: Mask::new(),
            calib: Cu::new(),
            src: Source::new(n_src as i32, pupil_size, 1024),
            atm_pssn: PSSn::new(),
            pssn: PSSn::new(),
            step: 0,
            atm: atm,
            d_mean_c: Cu::new(),
            x: Cu::new(),
        }
    }
    pub fn default(atm: &'a mut Atmosphere) -> Self {
        Self {
            sys: System::new(25.5, 4, 48, 16),
            lenslet_mask: Mask::new(),
            calib: Cu::new(),
            src: Source::new(1, 25.5, 1024),
            atm_pssn: PSSn::new(),
            pssn: PSSn::new(),
            step: 0,
            atm: atm,
            d_mean_c: Cu::new(),
            x: Cu::new(),
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
        println!("GLAO centroids #: {}", self.sys.wfs.n_centroids);

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
        println!("Centroid mask: [{}], nnz: {}", mask.len(), nnz);

        self.lenslet_mask.build(m);
        let mut f: ceo::Cu<f32> = ceo::Cu::vector(m);
        let mut q = mask.iter().cloned().map(|x| x as f32).collect::<Vec<f32>>();
        f.to_dev(&mut q);
        self.lenslet_mask.filter(&mut f);

        self
    }
    pub fn build_science_field(&mut self, src_zen: Vec<f32>, src_azi: Vec<f32>) -> &mut Self {
        let (_, gmt, _) = self.sys.devices();

        let n_src = src_zen.len();
        self.src.build("Vs", src_zen, src_azi, vec![0f32; n_src]);

        gmt.reset();
        self.src.through(gmt).xpupil();

        self.pssn.build(&mut self.src);
        self.atm_pssn.build(&mut self.src);

        self
    }
    pub fn build_single_layer_atmosphere(&mut self) -> &mut Self {
        self.atm.gmt_build(self.pssn.r0(), self.pssn.oscale);
        self
    }
    pub fn build_atmosphere(&mut self, altitude: Vec<f32>, xi0: Vec<f32>) -> &mut Self {
        let n_layer = altitude.len();
        self.atm.build(
            self.pssn.r0(),
            self.pssn.oscale,
            n_layer as i32,
            altitude,
            xi0,
            vec![0f32; n_layer],
            vec![0f32; n_layer],
        );
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
        println!("Calibration & Inversion in {}s", now.elapsed().as_secs());
        println!(
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
    pub fn wrap_up(&mut self) -> (Vec<f32>, Vec<f32>, f64, f64, f64) {
        self.atm_pssn.peek();
        self.pssn.peek();
        let mut fwhm = ceo::Fwhm::new();
        fwhm.build(&mut self.src);
        let atm_fwhm_x = ceo::Fwhm::atmosphere(
            self.pssn.wavelength as f64,
            self.pssn.r0() as f64,
            self.pssn.oscale as f64,
        )
        .to_arcsec();
        let atm_fwhm = fwhm
            .from_complex_otf(&self.atm_pssn.telescope_error_otf())
            .to_arcsec();
        let glao_fwhm = fwhm
            .from_complex_otf(&self.pssn.telescope_error_otf())
            .to_arcsec();
        (
            self.atm_pssn.estimates.clone(),
            self.pssn.estimates.clone(),
            atm_fwhm_x,
            atm_fwhm,
            glao_fwhm,
        )
    }
}

impl<'a> Iterator for GlaoSys<'a> {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        let n_kl = self.sys.m2_n_mode;
        let (gs, gmt, wfs) = self.sys.devices();

        let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);

        wfs.reset();
        gs.through(gmt).xpupil().through(self.atm).through(wfs);
        wfs.process()
            .fold_into(&mut self.d_mean_c, &mut self.lenslet_mask);

        self.src
            .through(gmt)
            .xpupil()
            .through(self.atm)
            .through(&mut self.atm_pssn);

        self.calib.qr_solve_as_ptr(&mut self.x, &mut self.d_mean_c);
        let h_x = self.x.from_dev();
        let mut k = 0;
        for s in 0..7 {
            for a in 1..n_kl {
                kl_coefs[s][a] -= h_x[k] as f64;
                k += 1;
            }
        }

        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);

        self.src
            .through(gmt)
            .xpupil()
            .through(self.atm)
            .through(&mut self.pssn);

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
