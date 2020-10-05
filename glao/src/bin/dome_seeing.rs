use ceo::{Cu, Gmt, Propagation, Source};
use cirrus;
use std::time::Instant;

pub struct DomeSeeing {
    region: String,
    bucket: String,
    folder: String,
    case: String,
    keys: Result<Vec<String>, String>,
    opd: Option<Vec<Vec<f32>>>,
    step: usize,
    buffer: Cu<f32>,
}
impl DomeSeeing {
    async fn new(region: &str, bucket: &str, folder: &str, case: &str) -> Self {
        DomeSeeing {
            region: region.to_owned(),
            bucket: bucket.to_owned(),
            folder: folder.to_owned(),
            case: case.to_owned(),
            keys: cirrus::list(region, bucket, &format! {"{}/{}",folder,case}, None).await,
            opd: Some(vec![]),
            step: 0,
            buffer: Cu::new(),
        }
    }
    async fn load_opd(&mut self) -> &mut Self {
        match &self.keys {
            Ok(keys) => match cirrus::load(&self.region, &self.bucket, &keys).await {
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
        self
    }
}
impl Propagation for DomeSeeing {
    /// Ray traces a `Source` through `Gmt`, ray tracing stops at the exit pupil
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        if self.opd.is_some() {
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

#[tokio::main]
async fn main() {
    let mut src = Source::new(1, 25.5, 512);
    src.build("V", vec![0.0], vec![0.0], vec![0.0]);
    let mut gmt = Gmt::new();
    gmt.build(1, None);
    src.through(&mut gmt).xpupil();
    println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);

    let mut ds = DomeSeeing::new(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020",
        "b2019_0z_0az_cd_12ms",
    )
    .await;

    let keys = ds.keys.as_ref().unwrap();
    let n_cfd_keys = keys.len();
    println!("CFD keys #: {} ; [{:?}]", n_cfd_keys, keys[..10].to_vec());

    ds.keys.as_mut().unwrap().truncate(10);
    let now = Instant::now();
    ds.load_opd().await;
    println!(
        "Downloaded {} files in {}s",
        ds.opd.as_ref().unwrap().len(),
        now.elapsed().as_secs()
    );

    let mut wfe_rms: Vec<f32> = vec![];
    src.through(&mut gmt).xpupil().through(&mut ds);
    wfe_rms.push(src.wfe_rms_10e(-9)[0]);
    while let Some(_) = ds.next() {
        src.through(&mut gmt).xpupil().through(&mut ds);
        wfe_rms.push(src.wfe_rms_10e(-9)[0]);
    }
    println!("{} steps in {}ms", ds.step, now.elapsed().as_millis());
    //println!("opd size: {}", ds.opd.unwrap().len());
    println!("Dome seeing WFE RMS: {:?}", &wfe_rms);
}
