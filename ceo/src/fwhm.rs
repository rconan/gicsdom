use super::Source;
use libm::j0;
use rayon::prelude::*;
use roots::find_root_brent;
use std::f64;

pub struct Fwhm {
    r: Vec<f64>,
    n_otf: usize,
    wavelength: f64,
    pub upper_bracket: f64,
}
impl Fwhm {
    pub fn new() -> Self {
        Fwhm {
            r: vec![],
            n_otf: 0,
            wavelength: 0f64,
            upper_bracket: 100f64,
        }
    }
    pub fn build(&mut self, src: &mut Source) -> &mut Self {
        self.wavelength = src.wavelength() as f64;
        let n = src.pupil_sampling as usize;
        let width = src.pupil_size;
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
        self.r = Vec::with_capacity(n_otf * n_otf);
        for i in 0..n_otf {
            let x = u[i] * d;
            for j in 0..n_otf {
                let y = u[j] * d;
                self.r.push(x.hypot(y));
            }
        }
        self.n_otf = n_otf;
        self
    }
    pub fn from_complex_otf(&mut self, otf: &[f32]) -> Vec<f64> {
        otf.par_chunks(self.n_otf * self.n_otf * 2)
            .map(|o| {
                let scalar_eq = |e: f64| -> f64 {
                    let mut s = 0f64;
                    for (_r, _o) in self.r.iter().zip(o.chunks(2)) {
                        let q = f64::consts::PI * e * _r;
                        let g = j0(q);
                        s += f64::from(_o[0]) * (0.5f64 - g) + f64::from(_o[1]) * g;
                    }
                    s
                };
                let root = find_root_brent(0f64, self.upper_bracket, &scalar_eq, &mut 1e-6f64);
                match root {
                    Ok(root) => self.wavelength * root,
                    Err(e) => {
                        println!("FWHM: {}", e);
                        f64::NAN
                    }
                }
            })
            .collect::<Vec<f64>>()
    }
    pub fn atmosphere(wavelength: f64, r0: f64, oscale: f64) -> f64 {
        0.9759 * (wavelength / r0) * (1.0 - 2.183 * (r0 / oscale).powf(0.356)).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{pssn::TelescopeError as TE, Atmosphere, Conversion, Gmt, PSSn};
    use std::time::Instant;

    #[test]
    fn fwhm_atmosphere_x() {
        let mut src = Source::new(1, 25.5, 1024);
        src.build("Vs", vec![0.0], vec![0.0], vec![0.0]);
        let mut gmt = Gmt::new();
        gmt.build(1, None);
        let mut pssn: PSSn<TE> = PSSn::new();
        src.through(&mut gmt);
        pssn.build(&mut src);
        let mut fwhm = Fwhm::new();
        fwhm.build(&mut src);
        let atm_fwhm_x = Fwhm::atmosphere(500e-9, pssn.r0() as f64, pssn.oscale as f64).to_arcsec();
        let atm_fwhm_n = fwhm.from_complex_otf(&pssn.atmosphere_otf());
        println!(
            "Atm. FWHM [arcsec]: {:.3}/{:.3}",
            atm_fwhm_x,
            atm_fwhm_n[0].to_arcsec()
        );
    }

    #[test]
    fn fwhm_atmosphere_n() {
        let mut src = Source::new(1, 25.5, 1024);
        src.build("Vs", vec![0.0], vec![0.0], vec![0.0]);
        let mut gmt = Gmt::new();
        gmt.build(1, None);
        let mut pssn: PSSn<TE> = PSSn::new();
        src.through(&mut gmt);
        pssn.build(&mut src);
        let mut fwhm = Fwhm::new();
        fwhm.build(&mut src);
        let mut atm = Atmosphere::new();
        atm.build(
            pssn.r0(),
            pssn.oscale,
            1,
            vec![0f32],
            vec![1f32],
            vec![0f32],
            vec![0f32],
        );
        let atm_fwhm_x0 =
            Fwhm::atmosphere(500e-9, pssn.r0() as f64, pssn.oscale as f64).to_arcsec();
        let atm_fwhm_x1 = fwhm.from_complex_otf(&pssn.atmosphere_otf());
        let mut k = 0;
        let now = Instant::now();
        loop {
            src.through(&mut gmt)
                .xpupil()
                .through(&mut atm)
                .through(&mut pssn);
            k += 1;
            if k == 10 {
                break;
            };
            atm.reset();
        }
        let atm_fwhm_n = fwhm.from_complex_otf(&pssn.telescope_error_otf());
        println!(
            "Atm. FWHM [arcsec]: {:.3}/{:.3}/{:.3} in {}s",
            atm_fwhm_x0,
            atm_fwhm_x1[0].to_arcsec(),
            atm_fwhm_n[0].to_arcsec(),
            now.elapsed().as_secs()
        );
    }
}
