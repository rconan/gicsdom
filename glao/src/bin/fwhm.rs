use ceo::Conversion;
use libm::j0;
use rayon::prelude::*;
use roots::find_root_brent;
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use std::f64;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use glao::system::Cn2;
//use std::time::Instant;


#[derive(Deserialize, Serialize)]
struct PSSnData {
    values: Vec<f32>,
    otf: Vec<f32>,
    #[serde(skip_deserializing)]
    fwhm: Vec<f64>,
}

fn main() {
    let n: usize = 512;
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

    cn2_profiles.par_iter().for_each(|cn2_prof| {

        let n_sample: usize = 1000;
        let filename = format!(
            "Results/glao_pssn_{:04}cn2_{:04}",
            cn2_prof.idx, n_sample
        );
        let data_path = Path::new(&filename);

        //println!("Loading data ...");
        //let now = Instant::now();
        let mut pssn_3128: PSSnData = {
            let file = File::open(data_path).unwrap();
            let reader = BufReader::with_capacity(1_000_000, file);
            pickle::from_reader(reader).unwrap()
        };
        //println!("... in {}s", now.elapsed().as_secs());

        let fwhm_def = |e: f64, o: &[f32]| -> f64 {
            let mut s = 0f64;
            for (_r, &_o) in r.iter().zip(o.iter().step_by(2)) {
                let q = f64::consts::PI * e * _r;
                let g = 1f64 - 2f64 * j0(q);
                s += f64::from(_o) * g;
            }
            s
        };

        //println!("Finding root ...");
        //let now = Instant::now();
        pssn_3128.fwhm = Vec::with_capacity(21);
        pssn_3128.fwhm = pssn_3128
            .otf
            .par_chunks(n_otf * n_otf * 2)
            .map(|o| {
                let fun = |e: f64| -> f64 { fwhm_def(e, o) };
                //println!("Root bracket: [{},{}],{}", fun(0f64), fun(20f64), j0(0f64));
                let root = find_root_brent(0f64, 20f64, &fun, &mut 1e-3f64).unwrap();
                (500e-9_f64 * root).to_arcsec()
            })
            .collect::<Vec<f64>>();
        //println!("... in {}s", now.elapsed().as_secs());
        println!("FWHM: {:?}", pssn_3128.fwhm);

        //println!("Saving data ...");
        //let now = Instant::now();
        let file = File::create(data_path).unwrap();
        let mut writer = BufWriter::with_capacity(1_000_000, file);
        pickle::to_writer(&mut writer, &pssn_3128, true).unwrap();
        //println!("... in {}s", now.elapsed().as_secs());
    });
}
