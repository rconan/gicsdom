use ceo::{Gmt, PSSn, Source};
use cirrus;
use cfd;
use log::LevelFilter;
use serde::{Deserialize, Serialize};
use simple_logger::SimpleLogger;
use std::{env, time::Instant};

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CfdCases {
    #[serde(rename = "baseline 2020")]
    baseline_2020: Vec<String>,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Results {
    time: Vec<f64>,
    wfe_rms: Vec<f32>,
    pssn: Vec<f32>,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .with_module_level("cirrus", LevelFilter::Info)
        .init()
        .unwrap();

    let job_idx = env::var("AWS_BATCH_JOB_ARRAY_INDEX")?
        .parse::<usize>()
        .expect("AWS_BATCH_JOB_ARRAY_INDEX parsing failed!");
    let duration = 2; /*env::var("N_SAMPLE")?
                           .parse::<usize>()
    
                       .expect("N_SAMPLE parsing failed!");*/
    let upload_results = false;

    let cfd_case = &cfd::get_cases()?[job_idx];
    println!("CFD CASE: {} with {} duration", cfd_case, duration);

    let mut src = Source::new(1, 25.5, 769);
    src.build("V", vec![0.0], vec![0.0], vec![0.0]);
    let mut gmt = Gmt::new();
    gmt.build(1, None);
    src.through(&mut gmt).xpupil();
    println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    let mut pssn: PSSn<ceo::pssn::TelescopeError> = PSSn::new();
    pssn.build(&mut src);

    let mut ds = cfd::DomeSeeing::new(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020",
        &cfd_case,
        duration,
        Some(40),
    );
    let now = Instant::now();
    ds.get_keys().await?.load_opd().await?;
    let keys = &ds.keys;
    println!(
        "CFD keys #: {} ; [{},...,{}]",
        ds.n_keys,
        keys.as_slice().first().unwrap(),
        keys.as_slice().last().unwrap()
    );
    println!(
        "Downloaded {} files in {}s",
        ds.opd.len(),
        now.elapsed().as_secs()
    );

    let mut t: Vec<f64> = vec![];
    let mut wfe_rms: Vec<f32> = vec![];
    let now = Instant::now();
    loop {
        src.through(&mut gmt)
            .xpupil()
            .through(&mut ds)
            .through(&mut pssn);
        t.push(ds.current_time);
        wfe_rms.push(src.wfe_rms_10e(-9)[0]);
        if ds.next().is_none() {
            break;
        };
    }
    println!("{} steps in {}s", ds.n_step, now.elapsed().as_secs());

    if upload_results {
        let results = Results {
            time: t,
            wfe_rms: wfe_rms,
            pssn: pssn.peek().estimates.clone(),
        };
        let key = format!("{}/{}/dome_seeing.pkl", "Baseline2020", cfd_case);
        cirrus::dump("us-west-2", "gmto.modeling", &key, &results).await?;
    } else {
        //println!("opd size: {}", ds.opd.unwrap().len());
        let w = t
            .clone()
            .into_iter()
            .zip(wfe_rms.clone().into_iter())
            .collect::<Vec<(f64, f32)>>();
        println!("Dome seeing WFE RMS: {:#?}", w);
        println!("PSSn: {:?}", pssn.peek().estimates);
    }
    Ok(())
}
