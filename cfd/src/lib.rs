use serde::{Deserialize, Serialize};
use serde_yaml;
use std::fs::File;

pub mod dome_seeing;

pub use self::dome_seeing::DomeSeeing;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CfdCases {
    #[serde(rename = "baseline 2020")]
    baseline_2020: Vec<String>,
}
pub fn get_cases() -> Result<Vec<String>,Box<dyn std::error::Error>> {
    let file = File::open("cfd/CFD_CASES.yaml")?;
    let cfd_cases_2020: CfdCases = serde_yaml::from_reader(file)?;
    Ok(cfd_cases_2020.baseline_2020)
}
