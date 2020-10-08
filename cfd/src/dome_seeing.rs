use ceo::{Cu, Propagation, Source};
use cirrus;

pub struct DomeSeeing {
    pub region: String,
    pub bucket: String,
    pub folder: String,
    pub case: String,
    pub keys: Vec<String>,
    pub n_keys: usize,
    first: usize,
    pub opd: Vec<Vec<f32>>,
    n_src: usize,
    step: usize,
    pub n_step: usize,
    buffer: Cu<f32>,
    rate: usize,
    time: Vec<f64>,
    pub current_time: f64,
}
impl DomeSeeing {
    pub fn new(
        region: &str,
        bucket: &str,
        folder: &str,
        case: &str,
        n_src: usize,
        rate: Option<usize>,
    ) -> Self {
        DomeSeeing {
            region: region.to_owned(),
            bucket: bucket.to_owned(),
            folder: folder.to_owned(),
            case: case.to_owned(),
            keys: vec![],
            n_keys: 0,
            first: 0,
            opd: vec![],
            n_src,
            step: 0,
            n_step: 0,
            buffer: Cu::new(),
            rate: rate.unwrap_or(1),
            time: vec![],
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
        self.time = b;
        Ok(self)
    }
    pub async fn load_opd(
        &mut self,
        n_last: Option<usize>,
    ) -> Result<&mut Self, Box<dyn std::error::Error>> {
        self.first = self.n_keys - n_last.unwrap_or(self.n_keys);
        let keys = &self.keys[self.first..];
        self.opd = cirrus::load::<Vec<f32>>(&self.region, &self.bucket, keys).await?;
        self.n_step = self.opd.len() * self.rate;
        self.buffer = Cu::vector(self.opd[0].len() * self.n_src);
        self.buffer.malloc();
        Ok(self)
    }
}
impl Propagation for DomeSeeing {
    /// Ray traces a `Source` through `Gmt`, ray tracing stops at the exit pupil
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        let idx = self.step / self.rate;
        self.current_time = self.time[self.first + idx];
        src.add_same(&mut self.buffer.to_dev(&mut self.opd[idx]));
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
        if self.step < self.n_step {
            Some(self.step)
        } else {
            None
        }
    }
}
