use ceo::{ceo, pssn::TelescopeError, set_gpu, Builder, FIELDDELAUNAY21, PSSN, SOURCE};
use cirrus;
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use serde_yaml as yaml;
use std::collections::BTreeMap;
use std::fs::File;

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
#[derive(Debug, Serialize, Deserialize, Default)]
struct Results {
    v_band: Band,
    h_band: Band,
}
#[derive(Debug, Serialize, Deserialize, Default)]
struct Band {
    pssn: Vec<Vec<f32>>,
    psu: f32,
}

fn rbm_to_pssn(data: Data) -> Band {
    let mut fem: FemRbm = data.data.fem;
    let n_sample = fem.m1_rbm.len();
    //println!("# sample: {}", n_sample);

    let mut gmt = ceo!(GMT, set_m1_n_mode = [27]);
    let _src_ = FIELDDELAUNAY21::new().set_band("H");
    let mut src = _src_.clone().build();
    let mut pssn = PSSN::<TelescopeError>::new().set_source(_src_).build();
    //println!("r0: {}m", pssn.r0());
    src.through(&mut gmt).xpupil().through(&mut pssn);
    pssn.reset();

    let sampling_time = 0.5e-3_f64;
    let peek_time = 1_f64;
    let peek_lag = (peek_time / sampling_time).round() as usize;
    let k0: usize = 10_000;
    let step: usize = 50;

    //let n_step = (n_sample - k0) / step;
    //println!("# step: {}", n_step);
    let mut pssn_sec: Vec<Vec<f32>> = vec![];
    let rbm = fem
        .m1_rbm
        .drain(k0..)
        .step_by(step)
        .zip(fem.m2_rbm.drain(k0..).step_by(step));
    for (k, (m1_rbm, m2_rbm)) in rbm.enumerate() {
        gmt.update42(Some(&m1_rbm), Some(&m2_rbm), None, None);
        src.through(&mut gmt).xpupil().through(&mut pssn);
        if (k * step % peek_lag) == 0 {
            pssn.peek();
            pssn_sec.push(pssn.estimates.clone());
        }
    }
    pssn.peek();
    pssn_sec.push(pssn.estimates.clone());
    let psu = pssn.spatial_uniformity();
    Band {
        pssn: pssn_sec,
        psu: psu,
    }
}

#[tokio::main]
async fn main() {
    let n_gpu = 1 as usize;

    let filename = "/home/ubuntu/CFD/CFD_CASES.yaml";
    let r = File::open(filename).unwrap();
    let mut cfd_cases: BTreeMap<String, Vec<String>> = yaml::from_reader(r).unwrap();
    let mut b2020_cases = cfd_cases.remove("baseline 2020").unwrap();
    b2020_cases.truncate(1);
    println!("CFD CASES: {:?}", b2020_cases);

    let mut handle = vec![];
    let cfd_cases = b2020_cases.clone();
    for (k, c) in cfd_cases.into_iter().enumerate() {
        handle.push(tokio::spawn(async move {
            set_gpu((k % n_gpu) as i32);
            let data: Data = {
                let key = format!(
                    "{}/{}/MT_FSM_IO_FSM_MountCtrl_TT7.rs.pkl",
                    "Baseline2020", c
                );
                let r = cirrus::reader("us-east-2", "gmto.starccm", &key)
                    .await
                    .unwrap();
                pickle::from_reader(r).expect("Failed")
            };
            rbm_to_pssn(data)
        }))
    }
    for (h,c) in handle.iter_mut().zip(b2020_cases.iter()) {
        let h_band_pssn = h.await.unwrap();
        let key = format!(
            "{}/{}/MT_FSM_IO_FSM_MountCtrl_TT7.pssn.pkl",
            "Baseline2020", c
        );
        let v_band_pssn: Band = {
            let r = cirrus::reader("us-east-2", "gmto.starccm", &key)
                .await
                .unwrap();
            pickle::from_reader(r).expect("Failed")
        };
        let results = Results {
            v_band: v_band_pssn,
            h_band: h_band_pssn,
        };

        cirrus::dump("us-east-2", "gmto.starccm", &key, &results)
            .await
            .unwrap();
    }
}
