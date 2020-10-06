use ceo::{Cu, Gmt, PSSn, Propagation, Source};
use cirrus;
use serde::{Deserialize, Serialize};
use serde_yaml;
use std::boxed::Box;
use std::{env, fs::File, time::Instant};

pub struct DomeSeeing {
    region: String,
    bucket: String,
    folder: String,
    case: String,
    keys: Result<Vec<String>, Box<dyn std::error::Error>>,
    n_keys: usize,
    first: usize,
    opd: Option<Vec<Vec<f32>>>,
    step: usize,
    buffer: Cu<f32>,
    time: Vec<f64>,
    current_time: f64,
}
impl DomeSeeing {
    async fn new(region: &str, bucket: &str, folder: &str, case: &str) -> Self {
        DomeSeeing {
            region: region.to_owned(),
            bucket: bucket.to_owned(),
            folder: folder.to_owned(),
            case: case.to_owned(),
            keys: Ok(vec![]),
            n_keys: 0,
            first: 0,
            opd: None,
            step: 0,
            buffer: Cu::new(),
            time: vec![],
            current_time: 0f64,
        }
    }
    async fn get_keys(&mut self) -> Result<&mut Self, Box<dyn std::error::Error>> {
        self.keys = cirrus::list(
            &self.region,
            &self.bucket,
            &format! {"{}/{}/OPDData_OPD_Data_",self.folder,self.case},
            None,
        )
        .await;
        self.n_keys = self.keys.as_ref().unwrap().len();
        self.time = self
            .keys
            .as_ref()
            .unwrap()
            .iter()
            .map(|x| x.split('/').last().unwrap()[17..29].parse::<f64>().unwrap())
            .collect::<_>();
        Ok(self)
    }
    async fn load_opd(
        &mut self,
        n_last: Option<usize>,
    ) -> Result<&mut Self, Box<dyn std::error::Error>> {
        self.first = self.n_keys - n_last.unwrap_or(self.n_keys);
        match &self.keys {
            Ok(keys) => match cirrus::load(&self.region, &self.bucket, &keys[self.first..]).await {
                Ok(opd) => {
                    self.buffer = Cu::vector(opd[0].len());
                    self.buffer.malloc();
                    self.opd = Some(opd);
                }
                Err(e) => {
                    println!("Error: {}", e);
                    self.opd = None;
                }
            },
            Err(e) => {
                println!("Error: {}", e);
                self.opd = None;
            }
        };
        Ok(self)
    }
}
impl Propagation for DomeSeeing {
    /// Ray traces a `Source` through `Gmt`, ray tracing stops at the exit pupil
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        if self.opd.is_some() {
            self.current_time = self.time[self.first + self.step];
            let opd = &mut self.opd.as_mut().unwrap();
            let data = &mut opd[self.step];
            src.add(&mut self.buffer.to_dev(data));
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}

impl Iterator for DomeSeeing {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        if self.step + 1 < self.opd.as_ref().unwrap_or(&vec![]).len() {
            self.step += 1;
            Some(self.step)
        } else {
            None
        }
    }
}

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
    let job_idx = env::var("AWS_BATCH_JOB_ARRAY_INDEX")?.parse::<usize>()?;
    let n_sample = env::var("N_SAMPLE")?.parse::<usize>()?;

    let file = File::open("CFD_CASES.yaml")?;
    let cfd_cases_2020: CfdCases = serde_yaml::from_reader(file)?;
    let cfd_case = &cfd_cases_2020.baseline_2020[job_idx];
    println!("CFD CASE: {}", cfd_case);

    let mut src = Source::new(1, 25.5, 769);
    src.build("V", vec![0.0], vec![0.0], vec![0.0]);
    let mut gmt = Gmt::new();
    gmt.build(1, None);
    src.through(&mut gmt).xpupil();
    println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);
    let mut pssn: PSSn<ceo::pssn::TelescopeError> = PSSn::new();
    pssn.build(&mut src);

    let mut ds = DomeSeeing::new("us-west-2", "gmto.modeling", "Baseline2020", &cfd_case).await;
    let now = Instant::now();
    ds.get_keys().await?.load_opd(Some(n_sample)).await?;
    let keys = ds.keys.as_ref().unwrap();
    println!(
        "CFD keys #: {} ; [{},...,{}]",
        ds.n_keys,
        keys.as_slice().first().unwrap(),
        keys.as_slice().last().unwrap()
    );
    println!(
        "Downloaded {} files in {}s",
        ds.opd.as_ref().unwrap().len(),
        now.elapsed().as_secs()
    );

    let mut t: Vec<f64> = vec![];
    let mut wfe_rms: Vec<f32> = vec![];
    src.through(&mut gmt).xpupil().through(&mut ds);
    t.push(ds.current_time);
    wfe_rms.push(src.wfe_rms_10e(-9)[0]);
    while let Some(_) = ds.next() {
        src.through(&mut gmt)
            .xpupil()
            .through(&mut ds)
            .through(&mut pssn);
        t.push(ds.current_time);
        wfe_rms.push(src.wfe_rms_10e(-9)[0]);
    }
    println!("{} steps in {}ms", ds.step+1, now.elapsed().as_millis());
    /*
    //println!("opd size: {}", ds.opd.unwrap().len());
    let w = t
        .clone()
        .into_iter()
        .zip(wfe_rms.clone().into_iter())
        .collect::<Vec<(f64, f32)>>();
    println!("Dome seeing WFE RMS: {:#?}", w);
    println!("PSSn: {:?}", pssn.peek().estimates);
     */

    let results = Results {
        time: t,
        wfe_rms: wfe_rms,
        pssn: pssn.peek().estimates.clone(),
    };
    let key = format!("{}/{}/dome_seeing.pkl", "Baseline2020", cfd_case);
    cirrus::dump("us-west-2", "gmto.modeling", &key, &results).await?;

    Ok(())
}
