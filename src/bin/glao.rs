use gicsdom::ceo;
use serde_pickle as pickle;
use std::collections::BTreeMap;
use std::f32;
use std::fs::File;
use std::time::Instant;

fn atmosphere_pssn() {
    let npx = 512;
    let n_src: usize = 6;
    let mut src = ceo::Source::new(n_src as i32, 25.5, npx);
    let zen: Vec<f32> = (0..n_src)
        .map(|x| ceo::Conversion::from_arcmin(x as f32))
        .collect();
    let azi = vec![0.0; n_src];
    let mag = vec![0.0; n_src];
    src.build("V", zen, azi, mag);

    let mut gmt = ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(None);
    src.through(&mut gmt).xpupil();
    println!("WFE RMS: {:?}nm", src.wfe_rms_10e(-9)[0]);

    let mut pssn: ceo::GPSSn<ceo::AtmosphereTelescopeError> = ceo::GPSSn::new();
    pssn.build(&mut src);
    println!("PSSn: {}", pssn.reset(&mut src));

    let mut atm = ceo::Atmosphere::new();
    atm.gmt_build(pssn.r0(), pssn.oscale);

    let n_sample = 100;
    let now = Instant::now();
    for _ in 0..n_sample {
        src.through(&mut gmt).xpupil().through(&mut atm);
        pssn.accumulate(&mut src);
        //println!("WFE RMS: {:?}nm",src.wfe_rms_10e(-9));
        //println!("PSSn: {}",pssn.peek(&mut src));
        atm.reset();
    }
    println!(
        "ET: {:.3}mn",
        now.elapsed().as_millis() as f64 * 1e-3 / 60f64
    );
    println!("PSSn: {}", pssn.reset(&mut src));
}

struct System {
    gmt: ceo::Gmt,
    wfs: ceo::GeometricShackHartmann,
    gs: ceo::Source,
}
impl System {
    fn new(pupil_size: f64, n_sensor: i32, n_lenslet: i32, n_px_lenslet: i32) -> Self {
        let d = pupil_size / n_lenslet as f64;
        Self {
            gmt: ceo::Gmt::new(),
            wfs: ceo::GeometricShackHartmann::new(n_sensor, n_lenslet, n_px_lenslet, d),
            gs: ceo::Source::default(),
        }
    }
    fn gmt_build(&mut self, m1_mode_type: &str, m1_n_mode: usize, m2_n_mode: usize) -> &mut Self {
        self.gmt
            .build_m1(m1_mode_type, m1_n_mode)
            .build_m2(Some(m2_n_mode));
        self
    }
    fn wfs_build(
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
    fn wfs_calibrate(&mut self, threshold: f64) -> &mut Self {
        self.gs.through(&mut self.gmt).xpupil();
        self.wfs.calibrate(&mut self.gs, threshold);
        self
    }
    fn through(&mut self) -> &mut Self {
        self.gs
            .through(&mut self.gmt)
            .xpupil()
            .through(&mut self.wfs);
        self
    }
    fn through_atmosphere(&mut self, atm: &mut ceo::Atmosphere) -> &mut Self {
        self.gs
            .through(&mut self.gmt)
            .through(atm)
            .through(&mut self.wfs);
        self
    }
    fn process(&mut self) {
        self.wfs.process();
    }
    fn m2_mode_calibrate(&mut self) -> Vec<f32> {
        let n_mode: usize = self.gmt.m2_n_mode;
        print!("M2 KL[{}] calibration ...", n_mode);
        let now = Instant::now();
        let mut kl_coefs = vec![vec![0f64; n_mode]; 7];
        let mut calib: Vec<f32> = vec![];
        for s in 0..7 {
            for a in 1..n_mode {
                let mut f = |x: f64| -> Vec<f32> {
                    kl_coefs[s][a] = x;
                    let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
                    self.gmt.set_m2_modes(&mut m);
                    self.through().process();
                    self.wfs.centroids.clone()
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
        println!(" in {}ms", now.elapsed().as_millis());
        println!(
            "M2 KL calibration size: {}X{}",
            self.wfs.n_centroids,
            calib.len() / self.wfs.n_centroids as usize
        );
        calib
    }
}

fn main() {
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
        .wfs_calibrate(wfs_intensity_threshold).through();

    let calib = on_axis_sys.m2_mode_calibrate();

    on_axis_sys.gmt.set_m2_modes_ij(1, 1, 1e-6);
    on_axis_sys.through().process();

    let mut a = ceo::CuFloat::new();
    a.malloc(calib.len()).into(&mut calib.clone());
    let mut d_s = ceo::CuFloat::new();
    d_s.malloc(on_axis_sys.wfs.n_centroids as usize).into(&mut on_axis_sys.wfs.centroids.clone());

    a.qr(on_axis_sys.wfs.n_centroids);
    let mut x = ceo::CuFloat::new();
    x.malloc(7*(N_KL-1));
    a.qr_solve(&mut x, &mut d_s);
    let mut h_x = x.from();
    h_x.truncate(20);
    println!("x: {:?}",h_x);

    //on_axis_sys.gmt.reset();
    //on_axis_sys.through();
    println!("WFE RMS: {:?}nm",on_axis_sys.gs.wfe_rms_10e(-9));
    /*
    let mut atm = ceo::Atmosphere::new();
    {
        let mut pssn: ceo::GPSSn<ceo::AtmosphereTelescopeError> = ceo::GPSSn::new();
        pssn.build(&mut on_axis_sys.gs);
        atm.gmt_build(pssn.r0(), pssn.oscale);
    }
    on_axis_sys.through_atmosphere(&mut atm).process();
    println!("WFE RMS: {:?}nm",on_axis_sys.gs.wfe_rms_10e(-9));
     */

    let mut d_s = ceo::CuFloat::new();
    d_s.malloc(on_axis_sys.wfs.n_centroids as usize).into(&mut on_axis_sys.wfs.centroids.clone());
    a.qr_solve(&mut x, &mut d_s);
    let mut h_x = x.from();
    let mut kl_coefs = vec![vec![0f64; N_KL]; 7];
    kl_coefs[1][1] = 1e-6;
    let mut k = 0;
    for s in 0..7 {
        for a in 1..N_KL {
            kl_coefs[s][a] -= h_x[k] as f64;
            k+=1;
        }
    }
    let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
    on_axis_sys.gmt.set_m2_modes(&mut m);
    on_axis_sys.through();
    println!("WFE RMS: {:?}nm",on_axis_sys.gs.wfe_rms_10e(-9));

    on_axis_sys.gmt.reset();
    on_axis_sys.through();
    println!("WFE RMS: {:?}nm",on_axis_sys.gs.wfe_rms_10e(-9));

    let mut atm = ceo::Atmosphere::new();
    {
        let mut pssn: ceo::GPSSn<ceo::AtmosphereTelescopeError> = ceo::GPSSn::new();
        pssn.build(&mut on_axis_sys.gs);
        atm.gmt_build(pssn.r0(), pssn.oscale);
    }
    on_axis_sys.through_atmosphere(&mut atm).process();
    println!("WFE RMS: {:?}nm",on_axis_sys.gs.wfe_rms_10e(-9));
    println!("S sum: {:?}",on_axis_sys.wfs.centroids.clone().iter().sum::<f32>());

    let mut d_s = ceo::CuFloat::new();
    d_s.malloc(on_axis_sys.wfs.n_centroids as usize).into(&mut on_axis_sys.wfs.centroids.clone());
    a.qr_solve(&mut x, &mut d_s);
    let mut h_x = x.from();
    let mut kl_coefs = vec![vec![0f64; N_KL]; 7];
    let mut k = 0;
    for s in 0..7 {
        for a in 1..N_KL {
            kl_coefs[s][a] -= h_x[k] as f64;
            k+=1;
        }
    }
    let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
    on_axis_sys.gmt.set_m2_modes(&mut m);
    on_axis_sys.through_atmosphere(&mut atm);
    println!("WFE RMS: {:?}nm",on_axis_sys.gs.wfe_rms_10e(-9));
    /*
    let mut data = BTreeMap::new();
    data.insert("calib".to_string(),calib);
    data.insert("s".to_string(),on_axis_sys.wfs.centroids.clone());

    let filename = "data.pkl";
    let mut file = File::create(filename).unwrap();
    pickle::to_writer(&mut file, &data, true).unwrap();
    */
}
