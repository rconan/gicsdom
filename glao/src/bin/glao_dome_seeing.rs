use ceo::Conversion;
//use cirrus;
use cfd;
use glao::glao_sys::{GlaoSys, ScienceField};
use log;
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
    atm_fwhm_x: f64,
    fwhm: Vec<f64>,
    pssn: Vec<f32>,
}

enum Loop {
    Closed,
    Open,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .with_module_level("cirrus", LevelFilter::Info)
        .with_module_level("glao", LevelFilter::Info)
        .init()
        .unwrap();

    let job_idx = env::var("AWS_BATCH_JOB_ARRAY_INDEX")
        .expect("AWS_BATCH_JOB_ARRAY_INDEX env var missing")
        .parse::<usize>()
        .expect("AWS_BATCH_JOB_ARRAY_INDEX parsing failed");
    let index_offset = env::var("N_INDEX_OFFSET").unwrap_or("N_0".to_owned())[2..]
        .parse::<usize>()
        .expect("N_INDEX_OFFSET parsing failed");

    let glao_loop = match env::var("LOOP").expect("LOOP:[OPEN,CLOSED] env var missing").as_str() {
        "OPEN" => Some(Loop::Open),
        "CLOSED" => Some(Loop::Closed),
        _ => None,
    }
    .expect("LOOP is either OPEN or CLOSED");

    let duration = env::var("DURATION")
        .expect("DURATION env var missing")
        .parse::<usize>()
        .expect("DURATION parsing failed!");

    let science_field = env::var("SCIENCE").expect("SCIENCE:[ONAXIS,DELAUNAY21] env var mission");

    let rate = 40;
    let upload_results = true;

    let cfd_case = &cfd::get_cases()?[job_idx + index_offset];
    println!("CFD CASE: {} with {} duration", cfd_case, duration);

    let n_px = 769;

    let secz = 1f32 / 30f32.to_radians().cos();
    let turb_cn2_height = [275, 425, 1250, 4000, 8000, 13000]
        .iter()
        .map(|x| *x as f32 * secz)
        .collect::<Vec<f32>>();
    let turb_cn2_xi0 = [0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751];
    let wind_speed = [5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187];
    let wind_direction = [0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462];
    let mut atm = ceo::Atmosphere::new();
    let mut science = match science_field.as_str() {
        "ONAXIS" => ScienceField::on_axis("Vs", n_px, None),
        "DELAUNAY21" => ScienceField::delaunay_21("Vs", n_px, None),
        _ => panic!("SCIENCE is either ONAXIS or DELAUNAY21")
    };
    println!("Building the science field ...");
    science.build();
    println!("Building the atmosphere ...");
    /*
    atm.build(
        science.pssn.r0(),
        science.pssn.oscale,
        6,
        turb_cn2_height,
        turb_cn2_xi0.to_vec(),
        wind_speed.to_vec(),
        wind_direction.to_vec(),
    );
    */
    atm.raytrace_build(
        science.pssn.r0(),
        science.pssn.oscale,
        6,
        turb_cn2_height,
        turb_cn2_xi0.to_vec(),
        wind_speed.to_vec(),
        wind_direction.to_vec(),
        25.5,
        n_px as i32,
        20f32.from_arcmin(),
        20f32,
        Some("glao_fiducial_atmosphere.bin"),
        Some(20),
    );
    println!("Building the GLAO system ...");
    let mut glao_4gs = GlaoSys::default(&mut atm, &mut science);
    let n_kl = 70;
    glao_4gs.build(6f32.from_arcmin(), n_kl, 0.5);

    let mut ds = cfd::DomeSeeing::new(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020",
        &cfd_case,
        duration,
        Some(rate),
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

    //let (gs, gmt, wfs) = glao_4gs.sys.devices();
    //gs.through(gmt).xpupil().through(&mut ds).through(wfs);

    let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
    println!(
        "Running {} loop ...",
        match glao_loop {
            Loop::Open => "open",
            Loop::Closed => "closed",
        }
    );
    let now = Instant::now();
    //let wfe_0 = glao_4gs.get_science(&mut ds);
    //let ps0 = glao_4gs.science.src.phase().clone();
    let mut t: Vec<f64> = vec![];
    let mut wfe: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = vec![];
    match glao_loop {
        Loop::Open => loop {
            t.push(ds.current_time);
            glao_4gs.atm.secs = ds.current_time.clone();
            wfe.push(glao_4gs.get_science(&mut ds));
            match ds.next() {
                Some(p) => {
                    if p == rate {
                        println!("Step #{}: reset PSSn!", p);
                        glao_4gs.science.pssn.reset();
                    }
                    if p%800==0 {
                        println!(
                            "{:9.3} {:5.0} {:?} {:?}",
                            ds.current_time,
                            wfe.last().unwrap().0[0],
                            wfe.last().unwrap().1,
                            wfe.last().unwrap().2
                        );
                    }
                }
                None => break,
            };
        },
        Loop::Closed => {
            println!("Calibrating the GLAO system ...");
            glao_4gs.calibration();
            loop {
                t.push(ds.current_time);
                glao_4gs.atm.secs = ds.current_time.clone();
                wfe.push(glao_4gs.get_science(&mut ds));
                glao_4gs.closed_loop(&mut ds, &mut kl_coefs, 0.5);
                /*
                */
                match ds.next() {
                    Some(p) => {
                        if p == rate {
                            println!("Step #{}: reset PSSn!", p);
                            glao_4gs.science.pssn.reset();
                        }
                        if p%800==0 {
                            println!(
                                "{:9.3} {:5.0} {:?} {:?}",
                                ds.current_time,
                                wfe.last().unwrap().0[0],
                                wfe.last().unwrap().1,
                                wfe.last().unwrap().2
                            );
                        }
                    }
                    None => break,
                };
            }
        }
    }

    println!("{} steps in {}s", ds.n_step, now.elapsed().as_secs());

    if upload_results {
        let (atm_fwhm_x, fwhm) = glao_4gs.science.glao_wrap_up();
        let results = Results {
            time: t,
            wfe,
            atm_fwhm_x,
            fwhm,
            pssn: glao_4gs.science.pssn.peek().estimates.clone(),
        };
        let key = format!(
            "{}/{}/glao_{}_dome_seeing.pkl",
            "Baseline2020",
            cfd_case,
            match glao_loop {
                Loop::Open => "open",
                Loop::Closed => "closed",
            }
        );
        println!("Uploading results to: s3://gmto.modeling/{}", key);
        cirrus::dump("us-west-2", "gmto.modeling", &key, &results).await?;
    } else {
        //println!("opd size: {}", ds.opd.unwrap().len());
        println!(
            "WFE: {:#?}",
            wfe.iter()
                .zip(t.iter())
                .take(5)
                .map(|x| format!("[{:5.3}: {:6.0}]", x.1, (x.0).0[0]))
                .collect::<Vec<String>>()
        );
        println!("PSSn: {:.?}", glao_4gs.science.pssn.peek().estimates);
    }

    Ok(())
}
