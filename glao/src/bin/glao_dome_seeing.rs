use ceo::Conversion;
//use cirrus;
use cfd;
use glao::glao_sys::{GlaoSys, ScienceField};
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
    wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)>,
    pssn: Vec<f32>,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .with_module_level("cirrus", LevelFilter::Info)
        .with_module_level("glao", LevelFilter::Info)
        .init()
        .unwrap();

    let job_idx = env::var("AWS_BATCH_JOB_ARRAY_INDEX")?
        .parse::<usize>()
        .expect("AWS_BATCH_JOB_ARRAY_INDEX parsing failed!");
    let n_sample = 2; /*env::var("N_SAMPLE")?
                                            .parse::<usize>()
                      .expect("N_SAMPLE parsing failed!");*/
    let rate = 40;
    let upload_results = true;

    let cfd_case = &cfd::get_cases()?[job_idx];
    println!("CFD CASE: {} with {} sample", cfd_case, n_sample);

    let n_px = 769;

    let mut atm = ceo::Atmosphere::new();
    //let mut science =
    //    ScienceField::delaunay_21("Vs", n_px, None);
    let mut science = ScienceField::on_axis("Vs", n_px, None);
    science.build();
    let mut glao_4gs = GlaoSys::default(&mut atm, &mut science);
    let n_kl = 70;
    glao_4gs.build(6f32.from_arcmin(), n_kl, 0.5).calibration();

    let mut ds = cfd::DomeSeeing::new(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020",
        &cfd_case,
        4,
        Some(rate),
    );
    let now = Instant::now();
    ds.get_keys().await?.load_opd(Some(n_sample)).await?;
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

    //let (gs, gmt, wfs) = glao_4gs.sys.devices();
    //gs.through(gmt).xpupil().through(&mut ds).through(wfs);

    let mut kl_coefs = vec![vec![0f64; n_kl]; 7];

    let now = Instant::now();
    //let wfe_0 = glao_4gs.get_science(&mut ds);
    //let ps0 = glao_4gs.science.src.phase().clone();
    let mut t: Vec<f64> = vec![];
    let mut wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = vec![];
    loop {
        t.push(ds.current_time);
        wfe.push(glao_4gs.get_science(&mut ds));
        glao_4gs.closed_loop(&mut ds, &mut kl_coefs, 0.5);
        //println!(" {:5.0} {:?} {:?}", wfe.0[0], wfe.1, wfe.2);
        if ds.next().is_none() {
            break;
        };
        match ds.next() {
            Some(p) => {
                if p == rate {
                    println!("Step #{}: reset PSSn!",p);
                    glao_4gs.science.pssn.reset();
                }
            }
            None => break,
        };
    }
    let pssn = &mut glao_4gs.science.pssn;
    /*
    let mut fwhm = ceo::Fwhm::new();
    fwhm.build(&mut glao_4gs.science.src);
    fwhm.upper_bracket = 2f64 / pssn.r0() as f64;
    let glao_fwhm = fwhm.from_complex_otf(&pssn.telescope_error_otf());
    */
    //let ps = glao_4gs.science.src.phase().clone();
    /*
    println!(
        "WFE RMS: [{:5.0},{:5.0}]nm ; PSSn: {:.5}",
        wfe_0.0[0],
        wfe.0[0],
        glao_4gs.science.pssn.peek().estimates[0]
    );
    */
    println!("{} steps in {}s", ds.n_step, now.elapsed().as_secs());

    if upload_results {
        let results = Results {
            time: t,
            wfe: wfe,
            pssn: pssn.peek().estimates.clone(),
        };
        let key = format!("{}/{}/glao_dome_seeing.pkl", "Baseline2020", cfd_case);
        cirrus::dump("us-west-2", "gmto.modeling", &key, &results).await?;
    } else {
        //println!("opd size: {}", ds.opd.unwrap().len());
        println!("PSSn: {:.?}", pssn.peek().estimates);
        /*
        let w = t
            .clone()
            .into_iter()
            .zip(wfe_rms.clone().into_iter())
            .collect::<Vec<(f64, f32)>>();
        println!("Dome seeing WFE RMS: {:#?}", w);
        println!("PSSn: {:?}", pssn.peek().estimates);
        */
    }

    Ok(())
}
