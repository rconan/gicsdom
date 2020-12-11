use ceo::{ceo, pssn::TelescopeError, set_gpu, Builder,PSSN};
use oqueue;
use rayon;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use serde_yaml as yaml;
use std::collections::BTreeMap;
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

const DATAPATH: &str = "/mnt/fsx/Baseline2020";

fn main() {
    let n_gpu = 8 as usize;
    let n_thread = 24;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_thread)
        .build()
        .unwrap();
    let oqueue = oqueue::Sequencer::stderr();

    let filename = "/home/ubuntu/CFD/CFD_CASES.yaml";
    let f = File::open(filename).unwrap();
    let mut cfd_cases: BTreeMap<String, Vec<String>> = yaml::from_reader(f).unwrap();
    let b2020_cases = cfd_cases.get_mut("baseline 2020").unwrap();
    //    b2020_cases.truncate(60);
    let cases = &b2020_cases[1..60].to_vec();

    pool.install(|| {
        cases.par_iter().for_each(|c| {
            let task = oqueue.begin();
            let thread_id = pool.current_thread_index().unwrap();
            //        println!("Thread ID#{}",thread_id);
            let gpu_id = (thread_id % n_gpu) as i32;
            write!(
                task,
                "CFD CASES: {} [thread #{:02}, gpu #{}]",
                c, thread_id, gpu_id
            );
            set_gpu(gpu_id);
            rbm_to_pssn(task, c)
        })
    });
}

fn rbm_to_pssn(task: oqueue::Task, cfd_case: &str) {
    let data: Data = {
        let key = String::from(format!(
            "{}/{}/MT_FSM_IO_FSM_MountCtrl_TT7.rs.pkl",
            DATAPATH, cfd_case
        ));
        write!(task, ", Loading {}", cfd_case.clone());
        let f = File::open(&key).expect(&format!("Failed to open {}", key));
        let r = BufReader::with_capacity(1_000_000, f);
        pickle::from_reader(r).expect("Failed")
    };

    let mut fem: FemRbm = data.data.fem;
    //let n_sample = fem.m1_rbm.len();
    //println!("# sample: {}", n_sample);

    let mut gmt = ceo!(GMT, set_m1_n_mode = [27]);
    let mut src = ceo!(FIELDDELAUNAY21, set_band = ["H"]);
    src.through(&mut gmt).xpupil();
    let mut pssn = PSSN::<TelescopeError>::new().set_source(&src).build();
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
    write!(task, ", Running ...");
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

    let h_band_pssn = Band {
        pssn: pssn_sec,
        psu: psu,
    };

    let key = format!(
        "{}/{}/MT_FSM_IO_FSM_MountCtrl_TT7.pssn.pkl",
        DATAPATH, cfd_case
    );
    writeln!(task, ", Saving to {}", cfd_case);
    let v_band_pssn: Band = {
        let r = File::open(&key).unwrap();
        pickle::from_reader(r).expect("Failed")
    };
    let results = Results {
        v_band: v_band_pssn,
        h_band: h_band_pssn,
    };
    let mut file = File::create(&key).unwrap();
    //println!("Saving {}...", filename);
    pickle::to_writer(&mut file, &results, true).unwrap();
}
