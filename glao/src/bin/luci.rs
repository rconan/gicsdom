use ceo;
use ceo::Conversion;
use csv;
use rayon;
use rayon::prelude::*;
use serde::Serialize;
use serde_pickle as pickle;
use std::collections::BTreeMap;
use std::f32;
use std::fs::File;
use std::time::Instant;

use glao::system::{atmosphere_pssn, Cn2, GlaoField, System};

#[derive(Debug, Serialize, Default)]
struct Results {
    function_name: String,
    args_in: BTreeMap<String, Vec<f32>>,
    args_out: BTreeMap<String, Vec<f32>>,
}

fn glao_pssn(
    n_sample: usize,
    r_not: f32,
    l_not: f32,
    n_layer: usize,
    altitude: Vec<f32>,
    xi0: Vec<f32>,
    filename: &str,
) {
    let glao_field_reader = File::open("glao_field.pkl").expect("File not found!");
    let glao_field: GlaoField =
        pickle::from_reader(glao_field_reader).expect("File loading failed!");
    let n_src = glao_field.zenith_arcmin.len();
    println!("GLAO field ({} samples):", n_src);
    println!(" * zenith : {:?}", glao_field.zenith_arcmin);
    println!(" * azimuth: {:?}", glao_field.azimuth_degree);
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

    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let n_kl = 170;

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
    println!("GLAO centroids #: {}", glao_sys.wfs.n_centroids);

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
            println!("lenslet mask: {}", s);
            for k in 0..v.len() {
                v[k] += if x[k] > 0f32 { 1usize } else { 0usize };
            }
        });
    let mask = v
        .iter()
        .map(|&x| if x > 0 { 1u8 } else { 0u8 })
        .collect::<Vec<u8>>();
    let nnz = mask.iter().cloned().map(|x| x as usize).sum::<usize>();
    println!("Centroid mask: [{}], nnz: {}", mask.len(), nnz);

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

    let mut src = ceo::Source::new(n_src as i32, pupil_size, 512);
    src.build("Vs", src_zen, src_azi, vec![0f32; n_src]);
    gmt.reset();
    src.through(gmt).xpupil();

    //let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
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

    let mut d_mean_c: ceo::Cu<f32> = ceo::Cu::vector(calib.n_row());
    let mut mean_c = vec![0f32; calib.n_row()];
    let mut x = ceo::Cu::<f32>::vector(calib.n_col());
    x.malloc();
    //let now = Instant::now();
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
    pssn.peek().xotf();

    let mut file = File::create(filename).unwrap();
    pickle::to_writer(&mut file, &pssn, true).unwrap();

    //    let mut file = File::create("glao_pssn.pkl").unwrap();
    //  pickle::to_writer(&mut file, &data, true).unwrap();
}

fn main() {
    let mut cn2_reader = csv::Reader::from_path("glao_cn2.csv").unwrap();
    let mut cn2_profiles: Vec<Cn2> = vec![];
    for result in cn2_reader.deserialize() {
        cn2_profiles.push(result.unwrap());
    }
    let secz = 1f32 / 30f32.to_radians().cos();

    let n_sample: usize = 1000;

    let n_gpu = 8 as usize;
    let n_thread = 32;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_thread)
        .build()
        .unwrap();

    pool.install(|| {
        cn2_profiles.par_iter().for_each(|cn2_prof| {
            let thread_id = pool.current_thread_index().unwrap();
            ceo::set_gpu((thread_id % n_gpu) as i32);

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

            let filename = format!("Results/glao_pssn_{:04}cn2_{:04}.pkl", cn2_prof.idx, n_sample);
            glao_pssn(
                n_sample,
                atm_r0 as f32,
                atm_oscale,
                7,
                turb_cn2_height,
                turb_cn2_xi0.to_vec(),
                &filename,
            );
        })
    });
}
