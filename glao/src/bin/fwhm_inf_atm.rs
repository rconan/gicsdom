use ceo;
use ceo::Conversion;
use glao::system::Cn2;
use indicatif::ParallelProgressIterator;
use libm::j0;
use rayon::prelude::*;
use roots::find_root_brent;
use serde_pickle as pickle;
use std::f64;
use std::fs::File;
use std::io::BufWriter;

fn main() {
    let mut gmt = ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(Some(0));
    let mut gs = ceo::Source::new(1, 25.5, 1024);
    gs.build("Vs", vec![0f32], vec![0f32], vec![0f32]);
    gs.through(&mut gmt).xpupil();
    let mut pssn: ceo::PSSn<ceo::pssn::TelescopeError> =
        ceo::PSSn::from_r0_and_outerscale(10e-2, 25f32);
    pssn.build(&mut gs);
    let r0 = pssn.r0() as f64;
    let seeing = (0.9759 * 500e-9 / r0).to_arcsec();
    let atm_fwhm = ((0.9759 * 500e-9 / r0)
        * (1.0 - 2.183 * (r0 / pssn.oscale as f64).powf(0.356)).sqrt())
    .to_arcsec();
    pssn.integrate(&mut gs);
    println!("WFE RMS: {}; PSSN: {}", gs.wfe_rms_10e(-9)[0], pssn.peek());

    let n: usize = 1024;
    let width = 25.5_f64;
    let n_otf = 2 * n - 1;
    let d = width / (n - 1) as f64;

    let h = (n_otf - 1) / 2 + 1;
    let mut u: Vec<f64> = vec![];
    for k in 0..h {
        u.push(k as f64);
    }
    for k in 1..h {
        u.push(k as f64 - h as f64);
    }
    let mut r: Vec<f64> = Vec::with_capacity(n_otf * n_otf);
    for i in 0..n_otf {
        let x = u[i] * d;
        for j in 0..n_otf {
            let y = u[j] * d;
            r.push(x.hypot(y));
        }
    }

    let atm_otf = pssn.buffer_otf();

    let fwhm_def = |e: f64| -> f64 {
        let mut s = 0f64;
        for (_r, _o) in r.iter().zip(atm_otf.chunks(2)) {
            let q = f64::consts::PI * e * _r;
            let g = 1f64 - 2f64 * j0(q);
            s += f64::from(_o[0]) * g;
        }
        s
    };

    let root = find_root_brent(0f64, 20f64, &fwhm_def, &mut 1e-3f64).unwrap();
    let fwhm = (500e-9_f64 * root).to_arcsec();
    println!("seeing: {:.3}arcsec", seeing);
    println!("~FWHM : {:.3}arcsec", atm_fwhm);
    println!("FWHM  : {:.3}arcsec", fwhm);

    let file = File::create(format!("atmosphere_fwhm_inf_{}px_otf.pkl", n)).unwrap();
    let mut writer = BufWriter::new(file);
    pickle::to_writer(&mut writer, &atm_otf, true).unwrap();
}

fn luci() {
    let n: usize = 1024;
    let width = 25.5_f64;
    let n_otf = 2 * n - 1;
    let d = width / (n - 1) as f64;

    let h = (n_otf - 1) / 2 + 1;
    let mut u: Vec<f64> = vec![];
    for k in 0..h {
        u.push(k as f64);
    }
    for k in 1..h {
        u.push(k as f64 - h as f64);
    }
    let mut r: Vec<f64> = Vec::with_capacity(n_otf * n_otf);
    for i in 0..n_otf {
        let x = u[i] * d;
        for j in 0..n_otf {
            let y = u[j] * d;
            r.push(x.hypot(y));
        }
    }

    let mut cn2_reader = csv::Reader::from_path("glao_cn2.csv").unwrap();
    let mut cn2_profiles: Vec<Cn2> = vec![];
    for result in cn2_reader.deserialize() {
        cn2_profiles.push(result.unwrap());
    }

    let fwhm = cn2_profiles
        .par_iter()
        .progress_count(50)
        .map(|cn2_prof| {
            let mut gmt = ceo::Gmt::new();
            gmt.build_m1("bending modes", 0).build_m2(Some(0));
            let mut gs = ceo::Source::new(1, 25.5, n as i32);
            gs.build("Vs", vec![0f32], vec![0f32], vec![0f32]);
            gs.through(&mut gmt).xpupil();
            let atm_r0 = 0.9759 * 500e-9 / cn2_prof.dimm.from_arcsec();
            let mut pssn: ceo::PSSn<ceo::pssn::TelescopeError> =
                ceo::PSSn::from_r0_and_outerscale(atm_r0 as f32, 25f32);
            pssn.build(&mut gs);
            let atm_otf = pssn.buffer_otf();

            let fwhm_def = |e: f64| -> f64 {
                let mut s = 0f64;
                for (_r, _o) in r.iter().zip(atm_otf.chunks(2)) {
                    let q = f64::consts::PI * e * _r;
                    let g = 1f64 - 2f64 * j0(q);
                    s += f64::from(_o[0]) * g;
                }
                s
            };

            let root = find_root_brent(0f64, 20f64, &fwhm_def, &mut 1e-6f64).unwrap();
            (500e-9_f64 * root).to_arcsec()
        })
        .collect::<Vec<f64>>();

    let file = File::create(format!("atmosphere_fwhm_{}px.pkl", n)).unwrap();
    let mut writer = BufWriter::new(file);
    pickle::to_writer(&mut writer, &fwhm, true).unwrap();
}
