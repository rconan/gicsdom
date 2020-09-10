use ceo;
use ceo::Conversion;
use libm::j0;
use roots::find_root_brent;
use serde_pickle as pickle;
use std::f64;
use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

fn main() {
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

    let mut gmt = ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(Some(0));
    let mut gs = ceo::Source::new(1, 25.5, n as i32);
    gs.build("Vs", vec![0f32], vec![0f32], vec![0f32]);
    gs.through(&mut gmt).xpupil();
    let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> =
        ceo::PSSn::from_r0_and_outerscale(10e-2, 25f32);
    pssn.build(&mut gs);

    let root = {
        let mut pssn: ceo::PSSn<ceo::pssn::TelescopeError> =
            ceo::PSSn::from_r0_and_outerscale(10e-2, 25f32);
        pssn.build(&mut gs);
        gs.through(&mut gmt).xpupil();
        pssn.integrate(&mut gs);
        println!("PSSn: {}", pssn.peek());
        let otf = pssn.buffer_otf();
        let fwhm_def = |e: f64| -> f64 {
            let mut s = 0f64;
            for (_r, _o) in r.iter().zip(otf.chunks(2)) {
                let q = f64::consts::PI * e * _r;
                let g = 1f64 - 2f64 * j0(q);
                s += f64::from(_o[0]) * g;
            }
            s
        };
        find_root_brent(0f64, 20f64, &fwhm_def, &mut 1e-6f64).unwrap()
    };
    let fwhm_inf = (500e-9_f64 * root).to_arcsec();

    let mut atm = ceo::Atmosphere::new();
    atm.build(
        pssn.r0(),
        pssn.oscale,
        1i32,
        vec![0f32],
        vec![1f32],
        vec![0f32],
        vec![0f32],
    );

    let now = Instant::now();
    let n_sample = 10000;
    let mut fwhm_vs_sample = vec![];
    for k in 0..n_sample {
        gs.through(&mut gmt).xpupil().through(&mut atm);
        pssn.integrate(&mut gs);
        atm.reset();

        let m = if k < 1000 { 100 } else { 1000 };
        if k > 0 && (k + 1) % m == 0 {
            println!(
                "OTF integration ({:}) in {}ms",
                k+1,
                now.elapsed().as_millis()
            );

            let otf = pssn.telescope_error_otf();

            let fwhm_def = |e: f64| -> f64 {
                let mut sx = 0f64;
                for (_r, _o) in r.iter().zip(otf.chunks(2)) {
                    let q = f64::consts::PI * e * _r;
                    let g = 1f64 - 2f64 * j0(q);
                    sx += f64::from(_o[0]) * g;
                }
                sx
            };

            let root = find_root_brent(0f64, 10f64, &fwhm_def, &mut 1e-6f64).unwrap();
            let fwhm_n = (500e-9_f64 * root).to_arcsec();

            println!(
                "Inf. FWHM: {:.3} ; {:5} FWHM: {:.3} [arcsec]",
                fwhm_inf,
                k + 1,
                fwhm_n
            );
            fwhm_vs_sample.push((k+1,fwhm_n));
        }
    }

    let otf = pssn.telescope_error_otf();
    let file = File::create(format!("atmosphere_fwhm_n_{}px_sample.pkl", n)).unwrap();
    let mut writer = BufWriter::new(file);
    pickle::to_writer(&mut writer, &otf, true).unwrap();
}
