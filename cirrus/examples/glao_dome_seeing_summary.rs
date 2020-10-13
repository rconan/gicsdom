use cirrus;
use log::LevelFilter;
use serde::{Deserialize, Serialize};
use serde_yaml;
use simple_logger::SimpleLogger;
use std::fs::File;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CfdCases {
    #[serde(rename = "baseline 2020")]
    baseline_2020: Vec<String>,
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
struct Results {
    time: Vec<f64>,
    wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>,
    atm_fwhm_x: f64,
    fwhm: Vec<f64>,
    pssn: Vec<f32>,
}

#[tokio::main]
async fn main() {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .init()
        .unwrap();
    let file =
        File::open("/home/rconan/Dropbox/Documents/GMT/CFD/Python/CFD/CFD_CASES.yaml").unwrap();
    let cfd_cases: CfdCases = serde_yaml::from_reader(file).unwrap();

    println!(
        " {:<24} {:^6} {:>6} {:>8} {:>6} {:>6} {:>6} {:>9} {:>8}",
        "CFD case", "#", "dt", "T", "max", "min", "mean", "PSSn", "FWHM"
    );
    for cfd_case in cfd_cases.baseline_2020.iter().take(32) {
        let key = format!("{}/{}/glao_closed_dome_seeing.pkl", "Baseline2020", cfd_case);
        match cirrus::list("us-west-2", "gmto.modeling", &key, None).await {
            Ok(keys) => {
                let data = {
                    let data_v = cirrus::load::<Results>("us-west-2", "gmto.modeling", &keys)
                        .await
                        .unwrap();
                    data_v[0].clone()
                };
                let wfe_rms = data.wfe.iter().map(|x| (x.0)[0]).collect::<Vec<f32>>();
                let dt = data.time[1] - data.time[0];
                let duration = data.time.last().unwrap() - data.time[0];
                let n_sample = data.time.len();
                let wfe_rms_max = wfe_rms.iter().cloned().fold(-f32::INFINITY, f32::max);
                let wfe_rms_min = wfe_rms.iter().cloned().fold(f32::INFINITY, f32::min);
                let wfe_rms_mean: f32 =
                    (wfe_rms.iter().map(|x| x * x).sum::<f32>() / n_sample as f32).sqrt();
                println!(
                    " {:<24} {:>6} {:>6.3} {:>8.3} {:>6.0} {:>6.0} {:>6.0} {:>9.5} {:>8.3}",
                    cfd_case,
                    n_sample,
                    dt,
                    duration,
                    wfe_rms_max,
                    wfe_rms_min,
                    wfe_rms_mean,
                    data.pssn[12],
                    data.fwhm[12]
                );
            }
            Err(_) => {
                println!(" {:<24}", cfd_case);
            }
        }
    }
}
