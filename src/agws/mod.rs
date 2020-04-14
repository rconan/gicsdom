pub mod probe {

    use crate::astrotools;
    use crate::ceo;
    use crate::ceo::calibrations::{Mirror, Segment};
    use crate::ceo::LensletArray;
    use ceo::Conversion;
    use crossbeam_channel::{Receiver, Sender};
    use ndarray::Array2;
    use std::fmt;

    #[derive(Clone)]
    pub struct GmtState {
        pub m1_rbm: Vec<Vec<f64>>,
        pub m1_mode: Vec<Vec<f64>>,
        pub m2_rbm: Vec<Vec<f64>>,
    }
    impl GmtState {
        pub fn new(m1_n_mode: Option<usize>) -> GmtState {
            GmtState {
                m1_rbm: vec![vec![0.; 6]; 7],
                m1_mode: vec![vec![0.; m1_n_mode.or(Some(1)).unwrap()];7],
                m2_rbm: vec![vec![0.; 6]; 7],
            }
        }
    }
    impl GmtState {
        pub fn fill_from_array(
            &mut self,
            mirrors: Vec<Mirror>,
            segments: Vec<Vec<Segment>>,
            data: &Array2<f32>,
        ) -> &mut Self {
            let mut l = 0;
            for (sid, segment) in segments.iter().enumerate() {
                for mirror in mirrors.iter() {
                    for rbm in segment.iter() {
                        for a in rbm.range() {
                            match mirror {
                                Mirror::M1 => {
                                    self.m1_rbm[sid][a] = data[[l, 0]] as f64;
                                    l += 1;
                                }
                                Mirror::M1MODES => {
                                    self.m1_mode[sid][a] = data[[l, 0]] as f64;
                                    l += 1;
                                }
                                Mirror::M2 => {
                                    self.m2_rbm[sid][a] = data[[l, 0]] as f64;
                                    l += 1;
                                }
                            }
                        }
                    }
                }
            }
            self
        }
        pub fn acc_from_array(
            &mut self,
            mirrors: Vec<Mirror>,
            segments: Vec<Vec<Segment>>,
            data: &Array2<f32>,
        ) -> &mut Self {
            let mut l = 0;
            for (sid, segment) in segments.iter().enumerate() {
                for mirror in mirrors.iter() {
                    for rbm in segment.iter() {
                        for a in rbm.range() {
                            match mirror {
                                Mirror::M1 => {
                                    self.m1_rbm[sid][a] += data[[l, 0]] as f64;
                                    l += 1;
                                }
                                Mirror::M1MODES => {
                                    self.m1_mode[sid][a] += data[[l, 0]] as f64;
                                    l += 1;
                                }
                                Mirror::M2 => {
                                    self.m2_rbm[sid][a] += data[[l, 0]] as f64;
                                    l += 1;
                                }
                            }
                        }
                    }
                }
            }
            self
        }
    }
    impl fmt::Display for GmtState {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(
                f,
                "M1 |Txyz[nm],Rxyz[mas]| {:>14} M2 |Txyz[nm],Rxyz[mas]| {:>12} M1 BM [10e-9]\n{}",
                "","",
                (0..7)
                    .into_iter()
                    .map(|k| {
                        format!(
                            "|{},{}| |{},{}| |{} ...|",
                            self.m1_rbm[k][..3]
                                .iter()
                                .map(|x| format!("{:>+5.0}", x * 1e9))
                                .collect::<Vec<String>>()
                                .join(","),
                            self.m1_rbm[k][3..]
                                .iter()
                                .map(|x| format!("{:>+5.0}", x.to_mas()))
                                .collect::<Vec<String>>()
                                .join(","),
                            self.m2_rbm[k][..3]
                                .iter()
                                .map(|x| format!("{:>+5.0}", x * 1e9))
                                .collect::<Vec<String>>()
                                .join(","),
                            self.m2_rbm[k][3..]
                                .iter()
                                .map(|x| format!("{:>+5.0}", x.to_mas()))
                                .collect::<Vec<String>>()
                                .join(","),
                            self.m1_mode[k].
                                iter()
                                .map(|x| format!("{:>+5.0}", x * 1e9))
                                .take(10)
                                .collect::<Vec<String>>()
                                .join(","),
                        )
                    })
                    .collect::<Vec<String>>()
                    .join("\n")
            )
        }
    }

    #[derive(Clone)]
    pub struct Message {
        pub obs: astrotools::Observation,
        pub state: GmtState,
        pub opd: Option<Vec<f32>>,
    }

    pub trait Sensor {
        fn build(&mut self, zenith: Vec<f32>, azimuth: Vec<f32>, magnitude: Vec<f32>);
        fn through(&mut self, opd: Option<&mut ceo::CuFloat>) -> &mut Self;
        fn guide_star(&mut self) -> &mut ceo::Source;
        fn gmt(&mut self) -> &mut ceo::Gmt;
        fn detector(&mut self) -> &mut ceo::Imaging;
        fn lenslet_array(&self) -> LensletArray;
        fn pixel_scale(&mut self) -> f64;
        fn noise(&self) -> ceo::imaging::NoiseDataSheet;
    }

    pub struct SH48 {
        pub gmt: ceo::Gmt,
        m1_n_mode: usize,
        pub sensor: ceo::Imaging,
        pub optics: LensletArray,
        pub n_px_lenslet: i32,
        pub n_px_framelet: i32,
        binning: i32,
        pub guide_star: ceo::Source,
        guide_star_band: String,
        pub noise: ceo::imaging::NoiseDataSheet,
    }
    impl SH48 {
        pub fn new() -> SH48 {
            SH48 {
                gmt: ceo::Gmt::new(),
                m1_n_mode: 17,
                sensor: ceo::Imaging::new(),
                optics: LensletArray {
                    n_side_lenslet: 48,
                    lenslet_size: 25.5 / 48.,
                },
                n_px_lenslet: 16,
                n_px_framelet: 8,
                binning: 3,
                guide_star: ceo::Source::empty(),
                guide_star_band: String::from("R+I"),
                noise: ceo::imaging::NoiseDataSheet {
                    rms_read_out_noise: 0.5,
                    n_background_photon: 0.0,
                    noise_factor: 2f64.sqrt(),
                },
            }
        }
    }
    impl Sensor for SH48 {
        fn build(&mut self, zenith: Vec<f32>, azimuth: Vec<f32>, magnitude: Vec<f32>) {
            self.gmt.build(self.m1_n_mode, None);
            let dft_osf = 2;
            let n_px_imagelet = self.n_px_framelet * self.binning;
            self.sensor.build(
                1,
                self.optics.n_side_lenslet,
                self.n_px_lenslet,
                dft_osf,
                n_px_imagelet,
                self.binning,
            );
            self.guide_star = ceo::Source::new(
                1,
                self.optics.n_side_lenslet as f64 * self.optics.lenslet_size,
                self.optics.n_side_lenslet * self.n_px_lenslet + 1,
            );
            self.guide_star
                .build(&self.guide_star_band, zenith, azimuth, magnitude);
        }
        fn through(&mut self, opd: Option<&mut ceo::CuFloat>) -> &mut Self {
            if opd.is_none() {
                self.guide_star
                    .through(&mut self.gmt)
                    .xpupil()
                    .through(&mut self.sensor);
            }
            if opd.is_some() {
                self.guide_star
                    .through(&mut self.gmt)
                    .xpupil()
                    .add(opd.unwrap())
                    .through(&mut self.sensor);
            }
            self
        }
        fn guide_star(&mut self) -> &mut ceo::Source {
            &mut self.guide_star
        }
        fn gmt(&mut self) -> &mut ceo::Gmt {
            &mut self.gmt
        }
        fn detector(&mut self) -> &mut ceo::Imaging {
            &mut self.sensor
        }
        fn lenslet_array(&self) -> LensletArray {
            self.optics.clone()
        }
        fn pixel_scale(&mut self) -> f64 {
            0.5 * (self.binning as f64) * self.guide_star.wavelength() / self.optics.lenslet_size
        }
        fn noise(&self) -> ceo::imaging::NoiseDataSheet {
            self.noise
        }
    }
    impl Drop for SH48 {
        fn drop(&mut self) {
            drop(&self.gmt);
            drop(&self.sensor);
            drop(&self.guide_star);
        }
    }

    pub enum ProbeSensors {
        SH24,
        SH48,
    }

    pub struct Probe<'a, S: Sensor> {
        pub coordinates: astrotools::SkyCoordinates,
        guide_star_magnitude: f64,
        pub sensor: &'a mut S,
        exposure: u32,
        pub detector_frame: Vec<f32>,
        pub sensor_data0: ceo::Centroiding,
        pub sensor_data: ceo::Centroiding,
        sampling_time: f64,
    }
    impl<'a, S> Probe<'a, S>
    where
        S: Sensor,
    {
        pub fn new(
            coordinates: astrotools::SkyCoordinates,
            guide_star_magnitude: f64,
            sensor: &'a mut S,
            exposure: u32,
        ) -> Probe<S> {
            Probe {
                coordinates,
                guide_star_magnitude,
                sensor,
                exposure,
                detector_frame: vec![],
                sensor_data0: ceo::Centroiding::new(),
                sensor_data: ceo::Centroiding::new(),
                sampling_time: 0.0,
            }
        }
        pub fn init(&mut self, observation: &astrotools::Observation) {
            let (z, a) = self.coordinates.local_polar(observation);
            //println!("({},{})", z.to_degrees() * 60., a.to_degrees());
            self.sensor.build(
                vec![z as f32],
                vec![a as f32],
                vec![self.guide_star_magnitude as f32],
            );
            let resolution = self.sensor.detector().resolution() as usize;
            self.detector_frame = vec![0f32; resolution * resolution];
            let data_units = Some(self.sensor.pixel_scale());
            //println!("pixel size: {:}", data_units.unwrap().to_degrees() * 3600.);
            let data_units = Some(1f64); // TODO: fix units issue
            self.sensor_data0.build(
                self.sensor.lenslet_array().n_side_lenslet as u32,
                data_units,
            );
            self.sensor_data.build(
                self.sensor.lenslet_array().n_side_lenslet as u32,
                data_units,
            );
            self.sensor.guide_star().set_fwhm(3.16);
            self.sampling_time = observation.sampling_time;
        }
        pub fn calibrate_sensor(&mut self, intensity_threshold: f64) -> u32 {
            self.sensor.gmt().reset();
            self.sensor.detector().reset();
            self.sensor.through(None);
            self.sensor_data0.process(self.sensor.detector(), None);
            let n_valid_lenslet = self
                .sensor_data0
                .set_valid_lenslets(Some(intensity_threshold), None);
            /*
            println!(
                "# valid lenslet: {} ({:.0}%)",
                n_valid_lenslet,
                100f64 * n_valid_lenslet as f64
                    / (self.sensor.lenslet_array().n_side_lenslet
                        * self.sensor.lenslet_array().n_side_lenslet) as f64
            );
            */
            self.sensor.detector().reset();
            n_valid_lenslet
        }
        pub fn update(
            &mut self,
            observation: &astrotools::Observation,
            state: &GmtState,
        ) -> (f64, f64) {
            let (z, a) = self.coordinates.local_polar(observation);
            //println!("({},{})", z.to_degrees() * 60., a.to_degrees());
            //self.sensor.guide_star().update(vec![z], vec![a]);
            self.sensor.gmt().update(
                Some(&state.m1_rbm),
                Some(&state.m2_rbm),
                Some(&state.m1_mode),
            );
            (z, a)
        }
        pub fn through(&mut self, phase: Option<&mut ceo::CuFloat>) -> Option<Vec<f32>> {
            self.sensor.through(phase);
            let n_frame = self.sensor.detector().n_frame();
            //println!("{}WFE RMS: {:.0}nm ; frame #{:0.4}",
            //         termion::cursor::Goto(1,9),
            //         self.sensor.guide_star().wfe_rms_10e(-9)[0],n_frame);
            //print!("Probe sending centroids ...");
            if n_frame == self.exposure {
                let exposure = self.sampling_time * self.exposure as f64;
                //print!("exposure: {}s", exposure);
                let noise = Some(self.sensor.noise());
                self.sensor.detector().readout(exposure, noise);
                self.sensor_data
                    .process(self.sensor.detector(), Some(&self.sensor_data0));
                self.sensor.detector().reset();
                Some(
                    self.sensor_data
                        .grab()
                        .valids(Some(&self.sensor_data0.valid_lenslets)),
                )
            } else {
                None
            }
        }
    }
}