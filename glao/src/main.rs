use ceo;
use rayon;
use rayon::prelude::*;
use std::f32;
use std::time::Instant;

fn atmosphere_pssn(n_sample: usize, zen: Vec<f32>, azi: Vec<f32>) -> Vec<f32> {
    /*
    PSSn: [1.2741,1.2718,1.2664,1.2646,1.2633,1.2522] for
     * zen = [0,1,...,5]
     * azi = 0
     * n_sample = 100
    */
    let npx = 512;
    let n_src: usize = zen.len();
    let mut src = ceo::Source::new(n_src as i32, 25.5, npx);
    let mag = vec![0.0; n_src];
    src.build("V", zen, azi, mag);

    let mut gmt = ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(None);
    src.through(&mut gmt).xpupil();
    //println!("WFE RMS: {:?}nm", src.wfe_rms_10e(-9)[0]);

    let mut pssn: ceo::GPSSn<ceo::AtmosphereTelescopeError> = ceo::GPSSn::new();
    pssn.build(&mut src);
    //println!("PSSn: {}", pssn.reset(&mut src));

    let mut atm = ceo::Atmosphere::new();
    atm.gmt_build(pssn.r0(), pssn.oscale);

    //let now = Instant::now();
    for _ in 0..n_sample {
        src.through(&mut gmt).xpupil().through(&mut atm);
        pssn.accumulate(&mut src);
        //println!("WFE RMS: {:?}nm",src.wfe_rms_10e(-9));
        //println!("PSSn: {}",pssn.peek(&mut src));
        atm.reset();
    }
    /*
    println!(
        "ET: {:.3}mn",
        now.elapsed().as_millis() as f64 * 1e-3 / 60f64
    );
    println!("PSSn: {}", pssn.reset(&mut src));
    */
    pssn.reset(&mut src);
    pssn.estimates.clone()
}

#[allow(dead_code)]
struct System {
    gmt: ceo::Gmt,
    wfs: ceo::GeometricShackHartmann,
    gs: ceo::Source,
}
#[allow(dead_code)]
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
    fn m2_mode_calibrate(&mut self) -> ceo::Cu<f32> {
        let n_mode: usize = self.gmt.m2_n_mode;
        print!("M2 KL[{}] calibration ...", n_mode);
        let now = Instant::now();
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

        println!(" in {}ms", now.elapsed().as_millis());
        println!(
            "M2 KL calibration size: {}X{}",
            self.wfs.n_centroids,
            calib.len() / self.wfs.n_centroids as usize
        );
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
    let h_x = x.from_dev();
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
        let mut pssn: ceo::GPSSn<ceo::AtmosphereTelescopeError> = ceo::GPSSn::new();
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

fn main() {
    let n_sample = 100;
    let n_src = 6;
    let zen: Vec<f32> = (0..n_src)
        .map(|x| ceo::Conversion::from_arcmin(x as f32))
        .collect();
    let azi = vec![0.0; n_src];
    let za: Vec<(f32, f32)> = zen
        .iter()
        .zip(azi.iter())
        .map(|x| (*x.0, *x.1))
        .collect::<Vec<(f32, f32)>>();

    let n_gpu = 1 as usize;
    let n_thread = n_src;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_thread)
        .build()
        .unwrap();

    let now = Instant::now();
    let results = pool.install(|| {
        za.par_iter()
            .map(|x| {
                let thread_id = pool.current_thread_index().unwrap();
                ceo::set_gpu((thread_id % n_gpu) as i32);
                let (z, a) = x;
                atmosphere_pssn(n_sample, vec![*z], vec![*a])
            })
            .collect::<Vec<_>>()
    });
    println!("#{} sample in {}s", n_sample, now.elapsed().as_secs());
    println!("Results: {:?}", results);
}
