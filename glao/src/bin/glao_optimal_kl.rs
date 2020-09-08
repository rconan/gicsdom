use indicatif::ParallelProgressIterator;
use ceo::Conversion;
use glao::system::System;
use rayon::prelude::*;
use serde_pickle as pickle;
use std::collections::BTreeMap;
use std::f32;
use std::fs::File;
use std::time::Instant;

#[allow(unused_variables)]
fn glao(n_kl: usize, n_sample: usize) -> (f32, f32) {
    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    //let n_kl = 90;

    let n_gs = 4;
    let gs_zen = (0..n_gs).map(|_| 6f32.from_arcmin()).collect::<Vec<f32>>();
    let gs_azi = (0..n_gs)
        .map(|x| (x as f32) * 2f32 * f32::consts::PI / n_gs as f32)
        .collect::<Vec<f32>>();
    let mut glao_sys = System::new(pupil_size, n_gs as i32, n_lenslet, n_px_lenslet);
    glao_sys
        .gmt_build("bending modes", 27, n_kl)
        .wfs_build("Vs", gs_zen, gs_azi, vec![0f32; n_gs])
        .wfs_calibrate(wfs_intensity_threshold)
        .through();
    //println!("GLAO centroids #: {}", glao_sys.wfs.n_centroids);

    let gs = &mut glao_sys.gs;
    let gmt = &mut glao_sys.gmt;
    let wfs = &mut glao_sys.wfs;

    // CALIBRATION
    let now = Instant::now();
    let mut v = vec![0usize; (n_lenslet * n_lenslet) as usize];
    wfs.lenset_mask()
        .from_dev()
        .chunks((n_lenslet * n_lenslet) as usize)
        .for_each(|x| {
            let s = x
                .iter()
                .map(|x| if *x > 0f32 { 1 } else { 0 })
                .sum::<usize>();
            //println!("lenslet mask: {}", s);
            for k in 0..v.len() {
                v[k] += if x[k] > 0f32 { 1usize } else { 0usize };
            }
        });
    let mask = v
        .iter()
        .map(|&x| if x > 0 { 1u8 } else { 0u8 })
        .collect::<Vec<u8>>();
    let nnz = mask.iter().cloned().map(|x| x as usize).sum::<usize>();
    //println!("Centroid mask: [{}], nnz: {}", mask.len(), nnz);

    let n = 7 * (n_kl - 1);
    let m = (n_lenslet * n_lenslet) as usize;
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
        .chunks(m as usize)
        .map(|x| {
            let mut r = vec![];
            for k in 0..m {
                if mask[k] > 0 {
                    r.push(x[k])
                }
            }
            r
        })
        .flatten()
        .collect::<Vec<f32>>();
        //println!("Reduced calibration: {}", reduced_calib.len());

        calib.to_dev(&mut reduced_calib);
    }
    //full_calib.qr();
    calib.qr();
    /*
    println!("Calibration & Inversion in {}s", now.elapsed().as_secs());
    println!(
        "Calibration matrix size [{}x{}]",
        calib.n_row(),
        calib.n_col()
    );
     */

    let n_src = 1;
    let mut src = ceo::Source::new(n_src as i32, pupil_size, 512);
    src.build("Vs", vec![0f32], vec![0f32], vec![0f32]);
    gmt.reset();
    src.through(gmt).xpupil();

    //let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    pssn.build(&mut src);

    let mut atm = ceo::Atmosphere::new();
//    atm.gmt_build(pssn.r0(), pssn.oscale);
    atm.build(
        pssn.r0(),
        pssn.oscale,
        1i32,
        vec![500f32],
        vec![1f32],
        vec![0f32],
        vec![0f32],
    );

    let mut d_mean_c: ceo::Cu<f32> = ceo::Cu::vector(calib.n_row());
    let mut mean_c = vec![0f32; calib.n_row()];
    let mut x = ceo::Cu::<f32>::vector(calib.n_col());
    x.malloc();

    let mut a_wfe_var = 0f32;

    let mut ps0: Vec<f32> = vec![];
    let now = Instant::now();
    for _k_sample in 0..n_sample {
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
                let mut r: Vec<f32> = Vec::with_capacity(nnz);
                for k in 0..m {
                    if mask[k] > 0 {
                        r.push(x[k])
                    }
                }
                r
            })
            .flatten()
            .collect::<Vec<f32>>()
            .chunks(calib.n_row())
            .for_each(|c| {
                for k in 0..mean_c.len() {
                    mean_c[k] += c[k] / n_gs as f32;
                }
            });

        //let wfe_rms_0 = src.through(gmt).xpupil().through(&mut atm).wfe_rms_10e(-9)[0];
        //ps0 = src.phase().clone();

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

        let wfe_rms = src.wfe_rms_10e(-9)[0];
        a_wfe_var += wfe_rms * wfe_rms;
        /*
        if k_sample % 100 == 0 {
            println!(
                "#{:6}: WFE RMS: {:5.0}/{:5.0}nm ; PSSn: {}",
                k_sample,
                wfe_rms_0,
                wfe_rms,
                pssn.peek()
            );
        }
        */
        atm.reset()
    }
    //    println!("{} sample in {}s", n_sample, now.elapsed().as_secs());
    let a_wfe_rms = (a_wfe_var / n_sample as f32).sqrt();
    let pssn_val = pssn.peek().estimates[0];
    /*println!(
        "#{:03} WFE RMS: {:5.0}nm ; PSSn: {}",
        n_kl, a_wfe_rms, pssn_val
    );*/
    return (a_wfe_rms, pssn_val);
}

fn main() {
    let now = Instant::now();

    let n_gpu = 8 as usize;
    let results = (20..201)
        .step_by(10)
        .map(|x| x)
        .collect::<Vec<usize>>()
        .par_iter()
        .progress_count(10)
        .map(|n_kl| {
            let thread_id = rayon::current_thread_index().unwrap();
            ceo::set_gpu((thread_id % n_gpu) as i32);
            glao(*n_kl, 100)
        })
        .collect::<Vec<(f32, f32)>>();
    println!("Results: {:?}", results);

    //let results = glao(70, 2);
    //println!("Elapsed time: {}s", now.elapsed().as_secs());

    let mut data = BTreeMap::new();
    data.insert("results".to_owned(), results);
    let mut file = File::create("glao_optimal_kl.pkl").unwrap();
    pickle::to_writer(&mut file, &data, true).unwrap();
}
