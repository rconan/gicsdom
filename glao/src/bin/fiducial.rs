use ceo;
use ceo::Conversion;
use glao::glao_sys::{GlaoSys, ScienceField};
use indicatif::{ProgressBar, ProgressStyle};
use rayon;
use rayon::prelude::*;
use serde_pickle as pickle;
use std::f32;
use std::fs::File;
use std::time::Instant;

use glao::system::{Cn2, GlaoField, System};

#[allow(dead_code)]
/*
fn glao_pssn_fiducial(
    n_sample: usize,
    src_zen: Vec<f32>,
    src_azi: Vec<f32>,
) -> (Vec<f32>, Vec<f32>, f64, f64, f64) {
    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let n_kl = 70;

    let n_gs = 4;
    let gs_zen = (0..n_gs).map(|_| 6f32.from_arcmin()).collect::<Vec<f32>>();
    let gs_azi = (0..n_gs)
        .map(|x| (x as f32) * 2f32 * f32::consts::PI / n_gs as f32)
        .collect::<Vec<f32>>();
    let mut glao_sys = System::new(pupil_size, n_gs as i32, n_lenslet, n_px_lenslet);
    glao_sys
        .gmt_build("S12", 2, n_kl)
        .wfs_build("Vs", gs_zen, gs_azi, vec![0f32; n_gs])
        .wfs_calibrate(wfs_intensity_threshold)
        .through();
    println!("GLAO centroids #: {}", glao_sys.wfs.n_centroids);

    let gs = &mut glao_sys.gs;
    let gmt = &mut glao_sys.gmt;
    let wfs = &mut glao_sys.wfs;

    // CALIBRATION
    let now = Instant::now();
    let m = (n_lenslet * n_lenslet) as usize;
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

    let mut lenslet_mask = ceo::Mask::new();
    lenslet_mask.build(m);
    let mut f: ceo::Cu<f32> = ceo::Cu::vector(m);
    let mut q = mask.iter().cloned().map(|x| x as f32).collect::<Vec<f32>>();
    f.to_dev(&mut q);
    lenslet_mask.filter(&mut f);

    //let n = 7 * (n_kl - 1);
    let mut calib: ceo::Cu<f32> = {
        let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
        on_axis_sys
            .gmt_build("bending modes", 27, n_kl)
            .wfs_build("Vs", vec![0f32], vec![0f32], vec![0f32])
            .wfs_calibrate(0f64)
            .through();
        on_axis_sys.m2_mode_calibrate_data(&mut lenslet_mask)
    };
    //full_calib.qr();
    calib.qr();
    println!("Calibration & Inversion in {}s", now.elapsed().as_secs());
    println!(
        "Calibration matrix size [{}x{}]",
        calib.n_row(),
        calib.n_col()
    );

    let n_src = src_zen.len();
    let mut src = ceo::Source::new(n_src as i32, pupil_size, 1024);
    src.build("Vs", src_zen, src_azi, vec![0f32; n_src]);
    gmt.reset();
    src.through(gmt).xpupil();

    //let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    pssn.build(&mut src);
    let mut atm_pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    atm_pssn.build(&mut src);

    let mut atm = ceo::Atmosphere::new();
    atm.gmt_build(pssn.r0(), pssn.oscale);

    let mut d_mean_c: ceo::Cu<f32> = ceo::Cu::vector(calib.n_row());
    d_mean_c.malloc();
    let mut x = ceo::Cu::<f32>::vector(calib.n_col());
    x.malloc();
    //gmt.set_m1_modes_ij(0, 1, 1f64);
    //(1..7).for_each(|i| gmt.set_m1_modes_ij(i, 0, 1f64));
    //let now = Instant::now();

    let pb = ProgressBar::new(n_sample as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7}"),
    );
    for _k_sample in 0..n_sample {
        pb.inc(1);
        let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);

        wfs.reset();
        gs.through(gmt).xpupil().through(&mut atm).through(wfs);
        wfs.process().fold_into(&mut d_mean_c, &mut lenslet_mask);

        src.through(gmt).xpupil().through(&mut atm);
        atm_pssn.integrate(&mut src);

        calib.qr_solve_as_ptr(&mut x, &mut d_mean_c);
        let h_x = x.from_dev();
        let mut k = 0;
        for s in 0..7 {
            for a in 1..n_kl {
                kl_coefs[s][a] -= h_x[k] as f64;
                k += 1;
            }
        }

        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);

        src.through(gmt).xpupil().through(&mut atm);
        pssn.integrate(&mut src);

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

        atm.reset()
    }
    pb.finish();
    //println!("{} sample in {}s", n_sample, now.elapsed().as_secs());
    //println!("PSSn: {}", pssn.peek());
    atm_pssn.peek();
    pssn.peek();
    let mut fwhm = ceo::Fwhm::new();
    fwhm.build(&mut src);
    let atm_fwhm_x =
        ceo::Fwhm::atmosphere(pssn.wavelength as f64, pssn.r0() as f64, pssn.oscale as f64)
            .to_arcsec();
    let atm_fwhm = fwhm
        .from_complex_otf(&atm_pssn.telescope_error_otf())
        .to_arcsec();
    let glao_fwhm = fwhm
        .from_complex_otf(&pssn.telescope_error_otf())
        .to_arcsec();
    (
        atm_pssn.estimates.clone(),
        pssn.estimates.clone(),
        atm_fwhm_x,
        atm_fwhm,
        glao_fwhm,
    )
}
 */
#[allow(dead_code)]
fn on_axis_glao() {
    /*
    let (atm_pssn, pssn, atm_fwhm_x, atm_fwhm, glao_fwhm) =
        glao_pssn_fiducial(1, vec![0f32], vec![0f32]);
    println!(
        "PSSN: {:?}/{:?}={} ; FWHM [arcsec]: {:.3}/{:.3}/{:.3}",
        pssn,
        atm_pssn,
        pssn[0] / atm_pssn[0],
        atm_fwhm_x,
        atm_fwhm,
        glao_fwhm
    );
    */
    let mut atm = ceo::Atmosphere::new();
    let mut science = ScienceField::on_axis("Vs", 1024, None);
    science.build();
    atm.gmt_build(science.pssn.r0(), science.pssn.oscale);
    let mut glao_4gs = GlaoSys::default(&mut atm, &mut science);
    glao_4gs.build(6f32.from_arcmin(), 70, 0.5).calibration();

    let n_sample = 1000;
    let pb = ProgressBar::new(n_sample as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7}"),
    );
    let (atm_pssn, pssn, atm_fwhm_x, atm_fwhm, glao_fwhm) = {
        for k in &mut glao_4gs {
            pb.inc(1);
            if k == n_sample {
                break;
            }
        }
        pb.finish();
        science.wrap_up().results()
    };
    println!(
        "PSSN: {:?}/{:?}={} ; FWHM [arcsec]: {:.3}/{:.3}/{:.3}",
        pssn,
        atm_pssn,
        pssn[0] / atm_pssn[0],
        atm_fwhm_x,
        atm_fwhm[0],
        glao_fwhm[0]
    );
}

#[allow(dead_code)]
fn main() {
    //}luci() {
    let mut cn2_reader = csv::Reader::from_path("glao_cn2.csv").unwrap();
    let mut cn2_profiles: Vec<Cn2> = vec![];
    for result in cn2_reader.deserialize() {
        cn2_profiles.push(result.unwrap());
    }
    let secz = 1f32 / 30f32.to_radians().cos();

    let cn2_prof = &cn2_profiles[0];
    let turb_cn2_height = [40f32, 125f32, 350f32, 1500f32, 4000f32, 8000f32, 16000f32]
        .iter()
        .map(|x| x * secz)
        .collect::<Vec<f32>>();
    let turb_cn2_xi0 = [
        cn2_prof.m40,
        cn2_prof.m125,
        cn2_prof.m350,
        cn2_prof.m1500,
        cn2_prof.m4000,
        cn2_prof.m8000,
        cn2_prof.m16000,
    ];

    let atm_r0 = 0.9759 * 500e-9 / cn2_prof.dimm.from_arcsec();
    let atm_oscale = 25f32;
    println!("Atmosphere: r0={} , L0={}", atm_r0, atm_oscale);

    let n_px = (4.0 * 25.5
        / ceo::PSSn::<ceo::pssn::TelescopeError>::r0_at_z(atm_r0 as f32, 30f32.to_radians()))
    .ceil()
    .min(1024f32) as usize;
    println!("Source pupil sampling: {}", n_px);

    let mut atm = ceo::Atmosphere::new();
    let mut science =
        ScienceField::delaunay_21("Vs", n_px, Some((atm_r0 as f32, atm_oscale as f32)));
    science.build();
    atm.build(
        science.pssn.r0(),
        science.pssn.oscale,
        7,
        turb_cn2_height.clone(),
        turb_cn2_xi0.to_vec().clone(),
        vec![0f32; 7],
        vec![0f32; 7],
    );
    let mut glao_4gs = GlaoSys::default(&mut atm, &mut science);
    glao_4gs.build(6f32.from_arcmin(), 70, 0.5).calibration();

    let (atm_pssn, pssn, atm_fwhm_x, atm_fwhm, glao_fwhm) = {
        let n_sample = 1000;
        let pb = ProgressBar::new(n_sample as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}"),
        );
        let mut pssn_previous = 0f32;
        let mut tol_count = 0;
        let mut tol = 0f32;
        for k in 0..n_sample {
            pb.inc(1);
            glao_4gs.next();
            if k % 25 == 0 {
                let pssn_current = glao_4gs.peek()[12];
                tol = (pssn_previous - pssn_current).abs();
                pb.set_message(&format!("Tol.: {:>7.4}",tol));
                if tol < 1e-4 {
                    tol_count += 1;
                } else {
                    tol_count = 0;
                }
                pssn_previous = pssn_current;
            }
            if tol_count == 2 {
                science.pssn_nsample_tol = Some((k,tol));
                break;
            }
        }
        pb.finish();
        science.wrap_up().dump("luci.pkl").unwrap().results()
    };
    println!(
        "PSSN: {:?}/{:?}={} ; FWHM [arcsec]: {:.3}/{:?}/{:?}",
        pssn,
        atm_pssn,
        pssn[0] / atm_pssn[0],
        atm_fwhm_x,
        atm_fwhm,
        glao_fwhm
    );
}

/*
fn main() {
    let n_gpu = 8 as usize;
    let n_thread = 24;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_thread)
        .build()
        .unwrap();

    let o = (0..8).map(|x| x as f32).collect::<Vec<f32>>();

    let fid_pssn = pool.install(|| {
        o.par_iter()
            .map(|oo| {
                let thread_id = pool.current_thread_index().unwrap();
                ceo::set_gpu((thread_id % n_gpu) as i32);
                let n_sample: usize = 100;
                //let src_zen: Vec<f32> = (0..6).map(|x| (x as f32).from_arcmin() ).collect::<Vec<f32>>();
                //let src_azi: Vec<f32> = vec![oo*f32::consts::PI/4f32;6];
                let src_zen: Vec<f32> = vec![6f32.from_arcmin()];
                let src_azi: Vec<f32> = vec![oo * f32::consts::PI / 4f32];

                glao_pssn_fiducial(n_sample, src_zen, src_azi)
            })
            .collect::<Vec<_>>()
    });
    println!("PSSn: {:?}", fid_pssn);
}
*/
