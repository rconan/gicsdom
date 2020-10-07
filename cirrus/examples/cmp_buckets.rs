use cirrus;
use serde::{Deserialize, Serialize};
use serde_yaml;
use simple_logger::SimpleLogger;
use std::fs::File;
use log::LevelFilter;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CfdCases {
    #[serde(rename = "baseline 2020")]
    baseline_2020: Vec<String>,
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

    let mut in_out_cfd_keys = (0_usize, 0_usize);
    println!(
        " {:<24} {:>12} {:>12} {:>8}",
        "CFD CASE", "STARCCM", "MODELING", "%"
    );
    for cfd_case in cfd_cases.baseline_2020.iter().take(60) {
        let prefix = format!("Baseline2020/{}/OPDData_OPD_Data_", cfd_case);
        let in_cfd_keys = cirrus::list("us-east-2", "gmto.starccm", &prefix, Some(".csv.gz"))
            .await
            .unwrap();
        let out_cfd_keys = cirrus::list("us-west-2", "gmto.modeling", &prefix, None)
            .await
            .unwrap();
        let n_in_cfd_keys = in_cfd_keys.len();
        let n_out_cfd_keys = out_cfd_keys.len();
        let pct_complete = 100f64 * (n_out_cfd_keys as f64) / (n_in_cfd_keys as f64);
        in_out_cfd_keys.0 += n_in_cfd_keys;
        in_out_cfd_keys.1 += n_out_cfd_keys;
        println!(
            " {:<24} {:>12} {:>12} {:>8.2}",
            cfd_case, n_in_cfd_keys, n_out_cfd_keys, pct_complete
        );
    }
    println!("{:-<60}", "");
    let pct_complete = 100f64 * (in_out_cfd_keys.1 as f64) / (in_out_cfd_keys.0 as f64);
    println!(
        " {:<24} {:>12} {:>12} {:>8.2}",
        "Total", in_out_cfd_keys.0, in_out_cfd_keys.1, pct_complete
    );
}
