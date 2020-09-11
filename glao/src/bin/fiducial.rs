use ceo;
use ceo::Conversion;
use indicatif::{ProgressBar,ProgressStyle};
use rayon;
use rayon::prelude::*;
use std::f32;
use std::time::Instant;

use glao::system::System;

fn glao_pssn_fiducial(
    n_sample: usize,
    src_zen: Vec<f32>,
    src_azi: Vec<f32>,
) -> (Vec<f32>, Vec<f32>) {
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

    let n = 7 * (n_kl - 1);
    let mut calib: ceo::Cu<f32> = ceo::Cu::array(2 * nnz, n);
    {
        let mut reduced_calib = {
            let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
            on_axis_sys
                .gmt_build("bending modes", 27, n_kl)
                .wfs_build("Vs", vec![0f32], vec![0f32], vec![0f32])
                .wfs_calibrate(0f64)
                .through();
            on_axis_sys.m2_mode_calibrate()
        }
        .from_dev()
        .chunks(m)
        .map(|x| {
            let (left, _): (Vec<f32>, Vec<u8>) = x
                .iter()
                .zip(mask.iter())
                .filter(|y| y.1 > &0u8)
                .collect::<Vec<(&f32, &u8)>>()
                .into_iter()
                .unzip();
            left.into_iter()
        })
        .flatten()
        .collect::<Vec<f32>>();
        println!("Reduced calibration: {}", reduced_calib.len());

        calib.to_dev(&mut reduced_calib);
    }
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
    //let mut atm_pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    //atm_pssn.build(&mut src);

    let mut atm = ceo::Atmosphere::new();
    atm.gmt_build(pssn.r0(), pssn.oscale);

    let mut d_mean_c: ceo::Cu<f32> = ceo::Cu::vector(calib.n_row());
    let mut mean_c = vec![0f32; calib.n_row()];
    let mut x = ceo::Cu::<f32>::vector(calib.n_col());
    x.malloc();
    gmt.set_m1_modes_ij(0, 1, 1f64);
    (1..7).for_each(|i| gmt.set_m1_modes_ij(i, 0, 1f64));
    //let now = Instant::now();
    let pb = ProgressBar::new(n_sample as u64);
    pb.set_style(ProgressStyle::default_bar()
                 .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7}"));
    let atm_pssn = vec![1.2731141f32];
    for _k_sample in 0..n_sample {
        pb.inc(1);
        let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
        let mut b = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut b);

        wfs.reset();
        gs.through(gmt).xpupil().through(&mut atm).through(wfs);
        wfs.process();

        for x in &mut mean_c {
            *x = 0f32;
        }
        wfs.centroids
            .from_dev()
            .chunks(m as usize)
            .map(|x| {
                let (left, _): (Vec<f32>, Vec<u8>) = x
                    .iter()
                    .zip(mask.iter())
                    .filter(|y| y.1 > &0u8)
                    .collect::<Vec<(&f32, &u8)>>()
                    .into_iter()
                    .unzip();
                left.into_iter()
            })
            .flatten()
            .collect::<Vec<f32>>()
            .chunks(calib.n_row())
            .for_each(|c| {
                for k in 0..mean_c.len() {
                    mean_c[k] += c[k] / n_gs as f32;
                }
            });

        //src.through(gmt).xpupil().through(&mut atm);
        //atm_pssn.integrate(&mut src);

        d_mean_c.to_dev(&mut mean_c);
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
    //println!("{} sample in {}s", n_sample, now.elapsed().as_secs());
    //println!("PSSn: {}", pssn.peek());
    //atm_pssn.peek();
    pssn.peek();
    (atm_pssn, pssn.estimates.clone())
}

fn main() {
    let (atm_pssn, pssn) = glao_pssn_fiducial(1000, vec![0f32], vec![0f32]);
    println!("PSSN: {:?}/{:?}={}", pssn, atm_pssn, pssn[0] / atm_pssn[0]);
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
