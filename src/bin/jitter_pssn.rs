use csv::{Reader, Writer};
use ensure::Present;
use gicsdom::astrotools;
use gicsdom::ceo::{set_gpu, Conversion, Gmt, PSSn, Source};
use indicatif::{ProgressBar, ProgressStyle};
use s3_sync::{Bucket, Object, Region, S3};
use serde::Deserialize;
use serde_pickle as pickle;
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

fn main() {
    let s3 = S3::new(Region::UsEast2);
    let bucket = Present(Bucket::from_name(String::from("gmto.starccm")));
    let key = String::from("MT_FSM_IO_FSM_MountCtrl_TT7.pkl");
    println!("Loading data from {}", &key);
    let object = Object::from_key(&bucket, key);
    let mut body = Vec::new();
    let now = Instant::now();
    s3.get_body(&Present(object))
        .expect("object body")
        .read_to_end(&mut body)
        .unwrap();
    let r = Cursor::new(body);
    let mut data: Data = pickle::from_reader(r).expect("Failed");
    println!("Data loaded in {}s", now.elapsed().as_secs());
    //println!("Data: {:#?}",data);

    //let filename = "/home/rconan/DATA/test.pkl";
    //let file = File::open(filename).unwrap();
    //println!("Loading {}...", filename);
    let mut fem: FemRbm = data.data.fem;//pickle::from_reader(file).expect("Failed");
    let n_sample = fem.m1_rbm.len();
    println!("# sample: {}", n_sample);
    /*
    println!("M1 RBM: #{} [{}]\n{:#?}", fem.m1_rbm.len(), fem.m1_rbm[0].len(), fem.m1_rbm.last().unwrap());
    let m1_srbm: Vec<Vec<f64>> = fem.m1_rbm.last().unwrap().chunks(7).map(|x| x.to_vec()).collect();
    println!("Chunked vector: {:?}",m1_srbm);
    println!("M2 RBM: #{} [{}]\n{:#?}", fem.m2_rbm.len(), fem.m2_rbm[0].len(), fem.m2_rbm.last().unwrap());
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
    println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);

    let mut pssn = PSSn::new();
    pssn.build(&mut src);
    pssn.reset(&mut src);

    let sampling_time = 0.5e-3_f64;
    let peek_time = 1_f64;
    let peek_lag = (peek_time / sampling_time).round() as usize;
    println!("Peek lag: {}", peek_lag);
    let k0: usize = 780_000;
    let step: usize = 50;

    let n_step = (n_sample - k0) / step;
    println!("# step: {}", n_step);
    let bar = ProgressBar::new(n_step as u64);
    bar.set_style(ProgressStyle::default_bar().template("{eta} {bar:40.cyan/blue} {msg}"));

    for (k, (m1_rbm, m2_rbm)) in fem
        .m1_rbm
        .drain(k0..)
        .step_by(step)
        .zip(fem.m2_rbm.drain(k0..).step_by(step))
        .enumerate()
    {
        bar.inc(1);
        gmt.update42(Some(&m1_rbm), Some(&m2_rbm), None);
        src.through(&mut gmt).xpupil();
        pssn.accumulate(&mut src);
        if (k * step % peek_lag) == 0 {
            pssn.peek(&mut src);
            //println!("PSSn[{}]: {}", k, pssn);
            bar.set_message(&format!("{}", pssn));
        }
    }
    bar.finish();
    let psu = pssn.spatial_uniformity();
    println!("PSSn: {} - PSU: {:.3}%", pssn, psu);
}
