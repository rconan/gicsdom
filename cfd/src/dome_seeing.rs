use ceo::{Cu, Propagation, Source};
use cirrus;

pub struct DomeSeeing {
    pub region: String,
    pub bucket: String,
    pub folder: String,
    pub case: String,
    pub keys: Vec<String>,
    pub n_keys: usize,
    pub opd: Vec<Vec<f32>>,
    step: usize,
    pub n_step: usize,
    buffer: Cu<f32>,
    duration: usize,
    rate: usize,
    time: Vec<f64>,
    sampling_time: f64,
    pub current_time: f64,
}
impl DomeSeeing {
    pub fn new(
        region: &str,
        bucket: &str,
        folder: &str,
        case: &str,
        duration: usize,
        rate: Option<usize>,
    ) -> Self {
        DomeSeeing {
            region: region.to_owned(),
            bucket: bucket.to_owned(),
            folder: folder.to_owned(),
            case: case.to_owned(),
            keys: vec![],
            n_keys: 0,
            opd: vec![],
            step: 0,
            n_step: duration * rate.unwrap_or(1),
            buffer: Cu::new(),
            duration,
            rate: rate.unwrap_or(1),
            time: vec![],
            sampling_time: 0f64,
            current_time: 0f64,
        }
    }
    pub async fn get_keys(&mut self) -> Result<&mut Self, Box<dyn std::error::Error>> {
        let keys = cirrus::list(
            &self.region,
            &self.bucket,
            &format! {"{}/{}/OPDData_OPD_Data_",self.folder,self.case},
            None,
        )
        .await?;
        self.n_keys = keys.len();
        let time = keys
            .iter()
            .map(|x| x.split('/').last().unwrap()[17..29].parse::<f64>().unwrap())
            .collect::<Vec<f64>>();
        let mut sorter = keys
            .into_iter()
            .zip(time.into_iter())
            .collect::<Vec<(String, f64)>>();
        sorter.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let (a, b): (Vec<_>, Vec<_>) = sorter.iter().cloned().unzip();
        self.keys = a;
        self.sampling_time = (b[1] - b[0]) / self.rate as f64;
        self.time = b;
        Ok(self)
    }
    pub async fn load_opd(&mut self) -> Result<&mut Self, Box<dyn std::error::Error>> {
        let first = self.n_keys - self.duration - 1;
        let keys = &self.keys[first..];
        self.opd = cirrus::load::<Vec<f32>>(&self.region, &self.bucket, keys).await?;
        self.buffer = Cu::vector(self.opd[0].len());
        self.buffer.malloc();
        Ok(self)
    }
}
impl Propagation for DomeSeeing {
    /// Ray traces a `Source` through `Gmt`, ray tracing stops at the exit pupil
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        let idx = self.step / self.rate;
        let k = self.step % self.rate;
        if k > 0 {
            let alpha = k as f64 / self.rate as f64;
            let mut opd = self.opd[idx]
                .iter()
                .zip(self.opd[idx + 1].iter())
                .map(|x| (1f64 - alpha) as f32 * x.0 + alpha as f32 * x.1)
                .collect::<Vec<f32>>();
            src.add_same(&mut self.buffer.to_dev(&mut opd));
        } else {
            src.add_same(&mut self.buffer.to_dev(&mut self.opd[idx]));
        };
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}

impl Iterator for DomeSeeing {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        self.step += 1;
        self.current_time = self.step as f64 * self.sampling_time;
        if self.step < self.n_step {
            Some(self.step)
        } else {
            None
        }
    }
}
