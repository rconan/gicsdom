use ceo;
use ceo::Conversion;
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use std::fs::File;
use std::time::Instant;

#[derive(Debug, Deserialize, Serialize, Default)]
pub struct GlaoField {
    pub zenith_arcmin: Vec<f32>,
    pub azimuth_degree: Vec<f32>,
}


#[derive(Serialize, Deserialize, Default, Debug)]
#[serde(rename_all = "UPPERCASE")]
pub struct Cn2 {
    pub idx: usize,
    pub date: String,
    pub utc: String,
    pub m40: f32,
    pub m125: f32,
    pub m350: f32,
    pub m1500: f32,
    pub m4000: f32,
    pub m8000: f32,
    pub m16000: f32,
    pub dimm: f64,
    pub wdir: f64,
    pub wspd: f64,
}

pub fn atmosphere_pssn(
    n_sample: usize,
    r_not: f32,
    l_not: f32,
    n_layer: usize,
    altitude: Vec<f32>,
    xi0: Vec<f32>,
    filename: &str,
) {
    /*
    PSSn: [1.2741,1.2718,1.2664,1.2646,1.2633,1.2522] for
     * zen = [0,1,...,5]
     * azi = 0
     * n_sample = 100
     */

    let glao_field_reader = File::open("glao_field.pkl").expect("File not found!");
    let glao_field: GlaoField =
        pickle::from_reader(glao_field_reader).expect("File loading failed!");
    let n_src = glao_field.zenith_arcmin.len();
    //println!("GLAO field ({} samples):", n_src);
    //println!(" * zenith : {:?}", glao_field.zenith_arcmin);
    //println!(" * azimuth: {:?}", glao_field.azimuth_degree);
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

    let npx = 512;
    let mut src = ceo::Source::new(n_src as i32, 25.5, npx);
    let mag = vec![0.0; n_src];
    src.build("Vs", src_zen, src_azi, mag);

    let mut gmt = ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(None);
    src.through(&mut gmt).xpupil();
    //println!("WFE RMS: {:?}nm", src.wfe_rms_10e(-9)[0]);

    let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> =
        ceo::PSSn::from_r0_and_outerscale(r_not, l_not);
    pssn.build(&mut src);

    let mut atm = ceo::Atmosphere::new();
    atm.build(
        pssn.r0(),
        pssn.oscale,
        n_layer as i32,
        altitude,
        xi0,
        vec![0f32; n_layer],
        vec![0f32; n_layer],
    );

    //let now = Instant::now();
    for _ in 0..n_sample {
        src.through(&mut gmt).xpupil().through(&mut atm);
        pssn.integrate(&mut src);
        /*
        if k % 100 == 0 {
            println!(
                "{:6}: WFE RMS: {:5.0}nm ; PSSn: {}",
                k,
                src.wfe_rms_10e(-9)[0],
                pssn.peek()
            );
        }
        */
        atm.reset();
    }
    /*
    println!(
        "ET: {:.3}mn",
        now.elapsed().as_millis() as f64 * 1e-3 / 60f64
    );
    println!("PSSn: {}", pssn.reset(&mut src));
    */
    pssn.peek().xotf();
    let mut file = File::create(filename).unwrap();
    pickle::to_writer(&mut file, &pssn, true).unwrap();
}

#[allow(dead_code)]
pub struct System {
    pub gmt: ceo::Gmt,
    pub wfs: ceo::GeometricShackHartmann,
    pub gs: ceo::Source,
}
#[allow(dead_code)]
impl System {
    pub fn new(pupil_size: f64, n_sensor: i32, n_lenslet: i32, n_px_lenslet: i32) -> Self {
        let d = pupil_size / n_lenslet as f64;
        Self {
            gmt: ceo::Gmt::new(),
            wfs: ceo::GeometricShackHartmann::new(n_sensor, n_lenslet, n_px_lenslet, d),
            gs: ceo::Source::default(),
        }
    }
    pub fn gmt_build(
        &mut self,
        m1_mode_type: &str,
        m1_n_mode: usize,
        m2_n_mode: usize,
    ) -> &mut Self {
        self.gmt
            .build_m1(m1_mode_type, m1_n_mode)
            .build_m2(Some(m2_n_mode));
        self
    }
    pub fn wfs_build(
        &mut self,
        band: &str,
        zenith: Vec<f32>,
        azimuth: Vec<f32>,
        magnitude: Vec<f32>,
    ) -> &mut Self {
        self.wfs.build();
        self.gs = self.wfs.new_guide_stars();
        self.gs.build(band, zenith, azimuth, magnitude);
        self
    }
    pub fn wfs_calibrate(&mut self, threshold: f64) -> &mut Self {
        self.gs.through(&mut self.gmt).xpupil();
        self.wfs.calibrate(&mut self.gs, threshold);
        self
    }
    pub fn through(&mut self) -> &mut Self {
        self.gs
            .through(&mut self.gmt)
            .xpupil()
            .through(&mut self.wfs);
        self
    }
    pub fn through_atmosphere(&mut self, atm: &mut ceo::Atmosphere) -> &mut Self {
        self.gs
            .through(&mut self.gmt)
            .through(atm)
            .through(&mut self.wfs);
        self
    }
    pub fn process(&mut self) {
        self.wfs.process();
    }
    pub fn m2_mode_calibrate(&mut self) -> ceo::Cu<f32> {
        let n_mode: usize = self.gmt.m2_n_mode;
        //print!("M2 KL[{}] calibration ...", n_mode);
        //let now = Instant::now();
        let m = self.wfs.n_centroids as usize;
        let n = 7 * (n_mode - 1);
        let mut kl_coefs = vec![vec![0f64; n_mode]; 7];
        let mut calib: Vec<f32> = vec![];
        for s in 0..7 {
            for a in 1..n_mode {
                let mut f = |x: f64| -> Vec<f32> {
                    kl_coefs[s][a] = x;
                    let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
                    self.gmt.set_m2_modes(&mut m);
                    self.wfs.reset();
                    self.through().process();
                    self.wfs.centroids.from_dev()
                };
                let cp = f(1e-6);
                let cm = f(-1e-6);

                calib.append(
                    &mut cp
                        .iter()
                        .zip(cm.iter())
                        .map(|x| 0.5e6 * (x.0 - x.1))
                        .collect::<Vec<f32>>(),
                );

                kl_coefs[s][a] = 0f64;
                let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
                self.gmt.set_m2_modes(&mut m);
            }
        }
        self.wfs.reset();
        self.gmt.reset();
        /*
        println!(" in {}ms", now.elapsed().as_millis());
        println!(
            "M2 KL calibration size: {}X{}",
            self.wfs.n_centroids,
            calib.len() / self.wfs.n_centroids as usize
        );
        */
        let mut a = ceo::Cu::<f32>::array(m, n);
        a.to_dev(&mut calib);
        a
    }
    pub fn m2_mode_calibrate_data(&mut self, lenslet_mask: &mut ceo::Mask) -> ceo::Cu<f32> {
        let n_mode: usize = self.gmt.m2_n_mode;
        //print!("M2 KL[{}] calibration ...", n_mode);
        //let now = Instant::now();
        let n = 7 * (n_mode - 1);
        let mut kl_coefs = vec![vec![0f64; n_mode]; 7];
        let mut calib: Vec<f32> = vec![];
        for s in 0..7 {
            for a in 1..n_mode {
                let mut f = |x: f64| -> Vec<f32> {
                    kl_coefs[s][a] = x;
                    let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
                    self.gmt.set_m2_modes(&mut m);
                    self.wfs.reset();
                    self.through().process();
                    self.wfs.filter(lenslet_mask).from_dev()
                };
                let cp = f(1e-6);
                let cm = f(-1e-6);

                calib.append(
                    &mut cp
                        .iter()
                        .zip(cm.iter())
                        .map(|x| 0.5e6 * (x.0 - x.1))
                        .collect::<Vec<f32>>(),
                );

                kl_coefs[s][a] = 0f64;
                let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
                self.gmt.set_m2_modes(&mut m);
            }
        }
        self.wfs.reset();
        self.gmt.reset();
        /*
        println!(" in {}ms", now.elapsed().as_millis());
        println!(
            "M2 KL calibration size: {}X{}",
            self.wfs.n_centroids,
            calib.len() / self.wfs.n_centroids as usize
        );
         */
        let m = calib.len()/n;
        let mut a = ceo::Cu::<f32>::array(m, n);
        a.to_dev(&mut calib);
        a
    }
}

#[test]
fn cmd_estimate() {
    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
    const N_KL: usize = 20;
    on_axis_sys
        .gmt_build("bending modes", 27, N_KL)
        .wfs_build("V", vec![0f32], vec![0f32], vec![0f32])
        .wfs_calibrate(wfs_intensity_threshold)
        .through();

    let mut calib = on_axis_sys.m2_mode_calibrate();
    //println!("calib sum: {:}",calib.iter().sum::<f32>());

    on_axis_sys.gmt.set_m2_modes_ij(1, 1, 1e-6);
    on_axis_sys.wfs.reset();
    on_axis_sys.through().process();

    let now = Instant::now();
    calib.qr();
    println!("QR factorization in {}ms", now.elapsed().as_millis());
    let now = Instant::now();
    let mut x = calib.qr_solve(&mut on_axis_sys.wfs.centroids);
    println!("QR solve in {}ms", now.elapsed().as_millis());
    let mut h_x = x.from_dev();
    h_x.truncate(20);
    println!("x: {:?}", h_x);
    assert_eq!((h_x[19] * 1e6f32).round(), 1f32);
}

#[test]
fn atmosphere_correction() {
    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
    const N_KL: usize = 20;
    on_axis_sys
        .gmt_build("bending modes", 27, N_KL)
        .wfs_build("V", vec![0f32], vec![0f32], vec![0f32])
        .wfs_calibrate(wfs_intensity_threshold)
        .through();

    let mut calib = on_axis_sys.m2_mode_calibrate();
    //println!("calib sum: {:}",calib.iter().sum::<f32>());

    on_axis_sys.gmt.set_m2_modes_ij(1, 1, 1e-6);
    on_axis_sys.wfs.reset();
    on_axis_sys.through().process();

    let now = Instant::now();
    calib.qr();
    println!("QR factorization in {}ms", now.elapsed().as_millis());

    on_axis_sys.gmt.reset();
    on_axis_sys.through();
    println!("WFE RMS: {:?}nm", on_axis_sys.gs.wfe_rms_10e(-9));

    let mut atm = ceo::Atmosphere::new();
    {
        let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
        pssn.build(&mut on_axis_sys.gs);
        atm.gmt_build(pssn.r0(), pssn.oscale);
    }
    on_axis_sys.wfs.reset();
    on_axis_sys.through_atmosphere(&mut atm).process();
    println!("WFE RMS: {:?}nm", on_axis_sys.gs.wfe_rms_10e(-9)[0].round());
    //    println!("S sum: {:?}",on_axis_sys.wfs.centroids.from_dev().iter().sum::<f32>());

    let mut x = calib.qr_solve(&mut on_axis_sys.wfs.centroids);
    let h_x = x.from_dev();
    let mut kl_coefs = vec![vec![0f64; N_KL]; 7];
    let mut k = 0;
    for s in 0..7 {
        for a in 1..N_KL {
            kl_coefs[s][a] -= h_x[k] as f64;
            k += 1;
        }
    }
    let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
    on_axis_sys.gmt.set_m2_modes(&mut m);
    on_axis_sys.through_atmosphere(&mut atm);
    let wfe_rms = on_axis_sys.gs.wfe_rms_10e(-9)[0].round();
    println!("WFE RMS: {:?}nm", wfe_rms);
    assert_eq!(wfe_rms, 658f32)
}
