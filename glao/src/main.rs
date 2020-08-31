mod system;

use ceo;
use ceo::Conversion;
use rayon;
use rayon::prelude::*;
use serde::Serialize;
use serde_pickle as pickle;
use std::collections::BTreeMap;
use std::f32;
use std::fs::File;
use std::time::Instant;

use system::{atmosphere_pssn, System};

#[derive(Debug, Serialize, Default)]
struct Results {
    function_name: String,
    args_in: BTreeMap<String, Vec<f32>>,
    args_out: BTreeMap<String, Vec<f32>>,
}

#[allow(dead_code)]
pub fn uncorrected_atmosphere_pssn() {
    let n_sample = 100;
    let n_src = 1;
    let zen: Vec<f32> = (0..n_src)
        .map(|x| ceo::Conversion::from_arcmin(x as f32))
        .collect();
    let azi = vec![0.0; n_src];
    let za: Vec<(f32, f32)> = zen
        .iter()
        .zip(azi.iter())
        .map(|x| (*x.0, *x.1))
        .collect::<Vec<(f32, f32)>>();

    let n_gpu = 8 as usize;
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

    let mut args_in = BTreeMap::new();
    args_in.insert("# sample".to_string(), vec![n_sample as f32]);
    args_in.insert("zenith".to_string(), zen);
    args_in.insert("azimuth".to_string(), azi);
    let mut args_out = BTreeMap::new();
    args_out.insert(
        "pssn".to_string(),
        results.into_iter().flatten().collect::<Vec<f32>>(),
    );

    let data = Results {
        function_name: "atmosphere_pssn".to_string(),
        args_in: args_in,
        args_out: args_out,
    };

    let mut file = File::create("atmosphere_pssn_100.pkl").unwrap();
    pickle::to_writer(&mut file, &data, true).unwrap();
}

#[allow(dead_code)]
fn optimal_kl() {
    ceo::set_gpu(1);
    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let mut a_wfe_rms: Vec<f32> = vec![];
    for n_kl in (20..201).step_by(10) {
        println!("N_KL: {}", n_kl);

        let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
        on_axis_sys
            .gmt_build("bending modes", 27, n_kl)
            .wfs_build("V", vec![0f32], vec![0f32], vec![0f32])
            .wfs_calibrate(wfs_intensity_threshold)
            .through();

        let mut calib = on_axis_sys.m2_mode_calibrate();
        calib.qr();

        let mut atm = ceo::Atmosphere::new();
        {
            let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
            pssn.build(&mut on_axis_sys.gs);
            atm.gmt_build(pssn.r0(), pssn.oscale);
        }

        let n_sample = 10;
        let mut a_wfe_var = 0f32;
        for _ in 0..n_sample {
            on_axis_sys.gmt.reset();
            on_axis_sys.wfs.reset();
            on_axis_sys.through_atmosphere(&mut atm).process();
            //let wfe_rms_0 = on_axis_sys.gs.wfe_rms_10e(-9)[0];

            let mut x = calib.qr_solve(&mut on_axis_sys.wfs.centroids);
            let h_x = x.from_dev();
            let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
            let mut k = 0;
            for s in 0..7 {
                for a in 1..n_kl {
                    kl_coefs[s][a] -= h_x[k] as f64;
                    k += 1;
                }
            }
            let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
            on_axis_sys.gmt.set_m2_modes(&mut m);
            on_axis_sys.through_atmosphere(&mut atm);
            let wfe_rms = on_axis_sys.gs.wfe_rms_10e(-9)[0];
            a_wfe_var += wfe_rms * wfe_rms;
            //println!("WFE RMS: {}/{}nm", wfe_rms_0, wfe_rms);
            atm.reset()
        }
        a_wfe_rms.push((a_wfe_var / n_sample as f32).sqrt());
    }
    println!("A. WFE RMS: {:?}", a_wfe_rms);
}

#[allow(dead_code)]
fn glao_pssn(n_sample: usize) {
    ceo::set_gpu(1);

    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let n_kl = 170;

    let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
    on_axis_sys
        .gmt_build("bending modes", 27, n_kl)
        .wfs_build("Vs", vec![0f32], vec![0f32], vec![0f32])
        .wfs_calibrate(0f64)
        .through();
    let now = Instant::now();
    let mut calib = on_axis_sys.m2_mode_calibrate();

    calib.qr();
    println!("Calibration & Inversion in {}s", now.elapsed().as_secs());
    println!(
        "Calibration matrix size [{}x{}]",
        calib.n_row(),
        calib.n_col()
    );

    let n_gs = 3;
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
    println!("GLAO centroids #: {}", glao_sys.wfs.n_centroids);

    let gs = &mut glao_sys.gs;
    let gmt = &mut glao_sys.gmt;
    let wfs = &mut glao_sys.wfs;

    let mut v = vec![0f32;(n_lenslet*n_lenslet) as usize];
    wfs.lenset_mask().from_dev().chunks((n_lenslet*n_lenslet) as usize).for_each( |x| {
        let s = x.iter().map(|x| if *x>0f32 { 1 } else { 0 }).sum::<usize>();
        println!("lenslet mask: {}",s);
        for k in 0..v.len() {
            v[k] += if x[k]>0f32 { 1f32 } else { 0f32 };
        }
    });

    let mut src = ceo::Source::new(1, pupil_size, 512);
    src.build("Vs", vec![0f32], vec![0f32], vec![0f32]);
    gmt.reset();
    src.through(gmt).xpupil();

    let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
    pssn.build(&mut src);

    let mut atm = ceo::Atmosphere::new();
    atm.gmt_build(pssn.r0(), pssn.oscale);

    //    let mut data = BTreeMap::new();

    let now = Instant::now();
    for k_sample in 0..n_sample {
        gmt.reset();
        wfs.reset();
        gs.reset();
        gs.through(&mut atm).through(gmt).opd2phase().through(wfs);
        wfs.process();
        //  data.insert("gs phase".to_string(), gs.phase().clone());

        let s = wfs.centroids.from_dev();
        let mut mean_c = vec![0f32; calib.n_row()];
        for c in s.as_slice().chunks(calib.n_row()) {
            for k in 0..mean_c.len() {
                mean_c[k] += c[k];
            }
        }
        let n = (n_lenslet*n_lenslet) as usize;
        for k in 0..n {
            if v[k]>0f32 {
                mean_c[k] /= v[k];
                mean_c[k+n] /= v[k];
            }
        }

        //      data.insert("c".to_string(), s);

        src.reset();
        let wfe_rms_0 = src
            .through(&mut atm)
            .through(gmt)
            .opd2phase()
            .wfe_rms_10e(-9)[0];
        //    data.insert("src phase".to_string(), src.phase().clone());

        let mut d_mean_c: ceo::Cu<f32> = ceo::Cu::vector(calib.n_row());
        d_mean_c.to_dev(&mut mean_c);

        let mut x = calib.qr_solve(&mut d_mean_c);
        let h_x = x.from_dev();
        let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
        let mut k = 0;
        for s in 0..7 {
            for a in 1..n_kl {
                kl_coefs[s][a] -= h_x[k] as f64;
                k += 1;
            }
        }

        let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
        gmt.set_m2_modes(&mut m);

        src.reset();
        src.through(&mut atm).through(gmt).opd2phase();
        pssn.accumulate(&mut src);
        let wfe_rms = src.wfe_rms_10e(-9)[0];
        if k_sample % 100 == 0 {
            println!(
                "#{:6}: WFE RMS: {:5.0}/{:5.0}nm ; PSSn: {}",
                k_sample,
                wfe_rms_0,
                wfe_rms,
                pssn.peek()
            );
        }

        atm.reset()
    }
    println!("{} sample in {}s", n_sample, now.elapsed().as_secs());
    println!("PSSn: {}", pssn.peek());

    //    let mut file = File::create("glao_pssn.pkl").unwrap();
    //  pickle::to_writer(&mut file, &data, true).unwrap();
}

#[allow(dead_code)]
fn glao_test(n_sample: usize) {
    ceo::set_gpu(1);

    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let n_kl = 170;

    /*
    let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
    on_axis_sys
        .gmt_build("bending modes", 27, n_kl)
        .wfs_build("Vs", vec![0f32], vec![0f32], vec![0f32])
        .wfs_calibrate(0f64)
        .through();
    let now = Instant::now();
    let mut calib = on_axis_sys.m2_mode_calibrate();

    calib.qr();
    println!("Calibration & Inversion in {}s", now.elapsed().as_secs());
    println!(
        "Calibration matrix size [{}x{}]",
        calib.n_row(),
        calib.n_col()
    );
     */

    let n_gs = 3;
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
    println!("GLAO centroids #: {}", glao_sys.wfs.n_centroids);

    //let gs = &mut glao_sys.gs;
    //let gmt = &mut glao_sys.gmt;
    let wfs = &mut glao_sys.wfs;

    //let mut lenslet_mask = wfs.lenset_mask().from_dev().chunks((n_lenslet*n_lenslet) as usize);
    let mut v = vec![0usize;(n_lenslet*n_lenslet) as usize];
    wfs.lenset_mask().from_dev().chunks((n_lenslet*n_lenslet) as usize).for_each( |x| {
        let s = x.iter().map(|x| if *x>0f32 { 1 } else { 0 }).sum::<usize>();
        println!("lenslet mask: {}",s);
        for k in 0..v.len() {
            v[k] += if x[k]>0f32 { 1 } else { 0 }
        }
    });

    let mut data = BTreeMap::new();
    data.insert("mask".to_owned(), v);
    let mut file = File::create("valid_lenslet.pkl").unwrap();
    pickle::to_writer(&mut file, &data, true).unwrap();

}

fn main() {
    //uncorrected_atmosphere_pssn();
    glao_pssn(100);
    //let onaxis_pssn = atmosphere_pssn(1, vec![0f32], vec![0f32]);
    //println!("On-axis PSSN: {:?}",onaxis_pssn);
}
