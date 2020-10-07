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
    wfe_rms: Vec<f32>,
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
        " {:<24} {:^6} {:>6} {:>8} {:>6} {:>6} {:>6} {:>9}",
        "CFD case", "#", "dt", "T", "max", "min", "mean", "PSSn"
    );
    for cfd_case in cfd_cases.baseline_2020.iter().take(60) {
        let key = format!("{}/{}/dome_seeing.pkl", "Baseline2020", cfd_case);
        let data: Results = {
            let data_v: Vec<Results> = cirrus::load("us-west-2", "gmto.modeling", &[key]).await.unwrap();
            data_v[0].clone()
        };
        let dt = data.time[1] - data.time[0];
        let duration = data.time.last().unwrap() - data.time[0];
        let n_sample = data.time.len();
        let wfe_rms_max = data.wfe_rms.iter().cloned().fold(-f32::INFINITY, f32::max);
        let wfe_rms_min = data.wfe_rms.iter().cloned().fold(f32::INFINITY, f32::min);
        let wfe_rms_mean: f32 = (data.wfe_rms.iter().map(|x| x*x).sum::<f32>()/n_sample as f32).sqrt();
        println!(
            " {:<24} {:>6} {:>6.2} {:>8.2} {:>6.0} {:>6.0} {:>6.0} {:>9.5}",
            cfd_case, n_sample, dt, duration, wfe_rms_max, wfe_rms_min, wfe_rms_mean, data.pssn[0]
        );
    }
}
