use csv::{Reader, Writer};
use ensure::Present;
use gicsdom::ceo::{set_gpu, Conversion, Gmt, PSSn, Source};
use indicatif::{ProgressBar, ProgressStyle};
use rayon;
use rayon::prelude::*;
use s3_sync::{Bucket, Object, ObjectBodyMeta, Region, S3};
use serde::{Deserialize, Serialize};
use serde_pickle as pickle;
use serde_yaml as yaml;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Cursor, Read};
use std::time::Instant;

#[derive(Debug, Deserialize, Default)]
struct ScienceField {
    zenith: f64,
    azimuth: f64,
}

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
    pssn: Vec<Vec<f32>>,
    psu: f32,
}

fn rbm_to_pssn(cfd_case: &str) {
    let s3 = S3::new(Region::UsEast2);
    let bucket = Present(Bucket::from_name(String::from("gmto.starccm")));
    let use_s3_repo = true;

    let data = match use_s3_repo {
        true => {
            let key = String::from(format!(
                "Baseline2020/{}/MT_FSM_IO_FSM_MountCtrl_TT7.rs.pkl",
                cfd_case
            ));
            //print!("Downloading data from {}", &key);
            let object = Object::from_key(&bucket, key);
            let mut body = Vec::new();
            let now = Instant::now();
            s3.get_body(&Present(object))
                .expect("object body")
                .read_to_end(&mut body)
                .unwrap();
            //println!(" in {}s", now.elapsed().as_secs());
            let r = Cursor::new(body);
            let now = Instant::now();
            //print!("Loading data ...");
            let data: Data = pickle::from_reader(r).expect("Failed");
            //println!(" in {}s", now.elapsed().as_secs());
            data
        }
        false => {
            let filename = "/home/ubuntu/MT_FSM_IO_FSM_MountCtrl_TT7.rs.pkl";
            let r = File::open(filename).unwrap();
            let now = Instant::now();
            //print!("Loading data ...");
            let data: Data = pickle::from_reader(r).expect("Failed");
            //println!(" in {}s", now.elapsed().as_secs());
            data
        }
    };

    let mut fem: FemRbm = data.data.fem; //pickle::from_reader(file).expect("Failed");
    let n_sample = fem.m1_rbm.len();
    //println!("# sample: {}", n_sample);
    /*
    //println!("M1 RBM: #{} [{}]\n{:#?}", fem.m1_rbm.len(), fem.m1_rbm[0].len(), fem.m1_rbm.last().unwrap());
    let m1_srbm: Vec<Vec<f64>> = fem.m1_rbm.last().unwrap().chunks(7).map(|x| x.to_vec()).collect();
    //println!("Chunked vector: {:?}",m1_srbm);
    //println!("M2 RBM: #{} [{}]\n{:#?}", fem.m2_rbm.len(), fem.m2_rbm[0].len(), fem.m2_rbm.last().unwrap());
    */

    let mut science_field: Vec<ScienceField> = vec![];
    let mut rdr = Reader::from_path("KPP_field_sampler.csv").unwrap();
    for result in rdr.deserialize() {
        science_field.push(result.unwrap());
    }
    let mut src = Source::new(science_field.len() as i32, 25.5, 512);
    let zen = science_field
        .iter()
        .map(|x| x.zenith.from_arcmin() as f32)
        .collect::<Vec<f32>>();
    let azi = science_field
        .iter()
        .map(|x| x.azimuth.to_radians() as f32)
        .collect::<Vec<f32>>();
    let mag = vec![0f32; science_field.len()];
    let star_color = "Vs";
    src.build(star_color, zen, azi, mag);

    let mut gmt = Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(None);
    src.through(&mut gmt).xpupil();
    //println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);

    let mut pssn = PSSn::new();
    pssn.build(&mut src);
    pssn.reset(&mut src);

    let sampling_time = 0.5e-3_f64;
    let peek_time = 1_f64;
    let peek_lag = (peek_time / sampling_time).round() as usize;
    //println!("Peek lag: {}", peek_lag);
    let k0: usize = 10_000;
    let step: usize = 50;

    let n_step = (n_sample - k0) / step;
/*
    //println!("# step: {}", n_step);
    let bar = ProgressBar::new(n_step as u64);
    bar.set_style(ProgressStyle::default_bar().template("{eta} {bar:40.cyan/blue} {msg}"));
     */
    let mut pssn_sec: Vec<Vec<f32>> = vec![];
    for (k, (m1_rbm, m2_rbm)) in fem
        .m1_rbm
        .drain(k0..)
        .step_by(step)
        .zip(fem.m2_rbm.drain(k0..).step_by(step))
        .enumerate()
    {
//        bar.inc(1);
        gmt.update42(Some(&m1_rbm), Some(&m2_rbm), None);
        src.through(&mut gmt).xpupil();
        pssn.accumulate(&mut src);
        if (k * step % peek_lag) == 0 {
            pssn.peek(&mut src);
            pssn_sec.push(pssn.estimates.clone());
            //println!("PSSn[{}]: {}", k, pssn);
            //bar.set_message(&format!("{}", pssn));
        }
    }
//    bar.finish();
    let psu = pssn.spatial_uniformity();
    pssn_sec.push(pssn.estimates.clone());
    //    println!("PSSn: {} - PSU: {:.3}%", pssn, psu);

    let results = Results {
        pssn: pssn_sec,
        psu: psu,
    };

    if use_s3_repo {
        let key = String::from(format!(
            "Baseline2020/{}/MT_FSM_IO_FSM_MountCtrl_TT7.pssn.pkl",
            cfd_case
        ));
        let data = pickle::to_vec(&results, true).unwrap();
        let body = Cursor::new(data);
        let object = Object::from_key(&bucket, key);
        s3.put_object(object, body, ObjectBodyMeta::default())
            .unwrap();
    } else {
        let filename = "/home/ubuntu/MT_FSM_IO_FSM_MountCtrl_TT7.pssn.pkl";
        let mut file = File::create(filename).unwrap();
        //println!("Saving {}...", filename);
        pickle::to_writer(&mut file, &results, true).unwrap();
    }
}

fn main() {
    let n_gpu = 8 as usize;
    let n_thread = 60;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_thread)
        .build()
        .unwrap();

    let filename = "/home/ubuntu/CFD/CFD_CASES.yaml";
    let r = File::open(filename).unwrap();
    let mut cfd_cases: BTreeMap<String, Vec<String>> = yaml::from_reader(r).unwrap();
    let b2020_cases = cfd_cases.get_mut("baseline 2020").unwrap();
    b2020_cases.truncate(60);
    println!("CFD CASES: {:?}", b2020_cases);

    pool.install(|| {
        b2020_cases.par_iter().for_each(|c| {
            let thread_id = pool.current_thread_index().unwrap();
            //        println!("Thread ID#{}",thread_id);
            set_gpu((thread_id % n_gpu) as i32);
            rbm_to_pssn(c)
        })
    });
}
