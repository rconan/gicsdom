use ceo::{ceo, pssn::TelescopeError, Builder, CEOType, PSSN, SOURCE, FIELDDELAUNAY21};
use cirrus;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use std::fs::File;
use std::io::BufReader;

#[derive(Debug, Deserialize, Default)]
struct FemRbm {
    #[serde(rename = "OSS_M1_lcl")]
    m1_rbm: Vec<Vec<f64>>,
    #[serde(rename = "MC_M2_lcl_6D")]
    m2_rbm: Vec<Vec<f64>>,
}
#[derive(Debug, Deserialize, Default)]
struct TimeSeries {
    #[serde(rename = "Time")]
    time: Vec<f64>,
    #[serde(rename = "Outputs")]
    fem: FemRbm,
}
#[derive(Debug, Deserialize, Default)]
struct Data {
    #[serde(rename = "SIMCEO")]
    simceo: String,
    #[serde(rename = "State space model")]
    state_space_model: String,
    #[serde(rename = "CFD case")]
    cfd_case: String,
    #[serde(rename = "Data")]
    data: TimeSeries,
}
#[derive(Debug, Serialize, Default)]
struct Results {
    v_band: Band,
    h_band: Band,
}
#[derive(Debug, Serialize, Default)]
struct Band {
    pssn: Vec<Vec<f32>>,
    psu: f32,
}

#[tokio::main]
async fn main() {
    let mut gmt = ceo!(GMT, set_m1_n_mode = [27]);
    let _src_ = FIELDDELAUNAY21::new();
    let mut src = _src_.clone().build();
    let mut pssn = PSSN::<TelescopeError>::new().set_source(_src_).build();
    src.through(&mut gmt).xpupil().through(&mut pssn);
    println!(
        "WFE RMS: {:?}nm ; PSSn: {}",
        src.wfe_rms_10e(-9),
        pssn.peek()
    );

    let use_s3 = true;
    let data: Data = if use_s3 {
        let cfd_case = "b2019_0z_0az_cd_12ms";
        let key = format!(
            "{}/{}/MT_FSM_IO_FSM_MountCtrl_TT7.rs.pkl",
            "Baseline2020", cfd_case
        );
        let r = cirrus::reader("us-east-2", "gmto.starccm", &key).await.unwrap();
        pickle::from_reader(r).expect("Failed")
    } else {
        let filename = "/home/ubuntu/DATA/MT_FSM_IO_FSM_MountCtrl_TT7.rs.pkl";
        let f = File::open(filename).unwrap();
        let r = BufReader::with_capacity(1000_0000, f);
        pickle::from_reader(r).expect("Failed")
    };
    let mut fem: FemRbm = data.data.fem;
    let n_sample = fem.m1_rbm.len();
    println!("# sample: {}", n_sample);

    let sampling_time = 0.5e-3_f64;
    let peek_time = 1_f64;
    let peek_lag = (peek_time / sampling_time).round() as usize;
    let k0: usize = 10_000;
    let step: usize = 50;

    let n_step = (n_sample - k0) / step;
    println!("# step: {}", n_step);
    let mut pssn_sec: Vec<Vec<f32>> = vec![];
    let rbm = fem
        .m1_rbm
        .drain(k0..)
        .step_by(step)
        .zip(fem.m2_rbm.drain(k0..).step_by(step));
    let bar = ProgressBar::new(n_step as u64);
    bar.set_style(ProgressStyle::default_bar().template("{eta} {bar:40.cyan/blue} {msg}"));
    for (k, (m1_rbm, m2_rbm)) in rbm.enumerate() {
        bar.inc(1);
        gmt.update42(Some(&m1_rbm), Some(&m2_rbm), None, None);
        src.through(&mut gmt).xpupil();
        pssn.accumulate(&mut src);
        if (k * step % peek_lag) == 0 {
            pssn.peek();
            pssn_sec.push(pssn.estimates.clone());
            //println!("PSSn[{}]: {}", k, pssn);
            bar.set_message(&format!("{}", pssn.estimates[0]));
        }
    }
    bar.finish();
    let psu = pssn.spatial_uniformity();
    pssn_sec.push(pssn.estimates.clone());
    let results = Results {
        v_band: Band {
            pssn: pssn_sec.clone(),
            psu: psu,
        },
        h_band: Band {
            pssn: pssn_sec,
            psu: psu,
        },
    };
    let filename = "jitter.pssn.pkl";
    let mut file = File::create(filename).unwrap();
    pickle::to_writer(&mut file, &results, true).unwrap();
}
