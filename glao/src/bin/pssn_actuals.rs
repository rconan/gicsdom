use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use std::f64;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use glao::system::Cn2;

#[derive(Deserialize, Serialize)]
struct PSSnData {
    values: Vec<f32>,
    otf: Vec<f32>,
    fwhm: Vec<f64>,
    #[serde(skip_deserializing)]
    actuals: Vec<f64>,
}

fn main() {
    let mut cn2_reader = csv::Reader::from_path("glao_cn2.csv").unwrap();
    let mut cn2_profiles: Vec<Cn2> = vec![];
    for result in cn2_reader.deserialize() {
        cn2_profiles.push(result.unwrap());
    }

    cn2_profiles.par_iter().progress_count(50).for_each(|cn2_prof| {
        let n_sample: usize = 1000;

        //println!("Loading data ...");
        //let now = Instant::now();
        let atmo_pssn: PSSnData = {
            let filename = format!(
                "Results/70KL/atmosphere_pssn_{:04}cn2_{:04}.pkl",
                cn2_prof.idx, n_sample
            );
            let data_path = Path::new(&filename);
            let file = File::open(data_path).unwrap();
            let reader = BufReader::with_capacity(1_000_000, file);
            pickle::from_reader(reader).unwrap()
        };

        let filename = format!("Results/70KL/glao_pssn_{:04}cn2_{:04}.pkl", cn2_prof.idx, n_sample);
        let data_path = Path::new(&filename);
        let mut glao_pssn: PSSnData = {
            let file = File::open(data_path).unwrap();
            let reader = BufReader::with_capacity(1_000_000, file);
            pickle::from_reader(reader).unwrap()
        };

        glao_pssn.actuals = glao_pssn
            .values
            .iter()
            .zip(atmo_pssn.values.iter())
            .map(|x| f64::from(*x.0) / f64::from(*x.1))
            .collect::<Vec<f64>>();

        let file = File::create(data_path).unwrap();
        let mut writer = BufWriter::with_capacity(1_000_000, file);
        pickle::to_writer(&mut writer, &glao_pssn, true).unwrap();
    });
}
