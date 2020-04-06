pub mod probe {

    use crate::astrotools;
    use crate::ceo;
    use crossbeam_channel::{bounded, unbounded, Receiver, Sender};

    pub trait Sensor {
        fn build(&mut self, zenith: Vec<f32>, azimuth: Vec<f32>, magnitude: Vec<f32>);
        fn through(&mut self) -> &mut Self;
        fn guide_star(&mut self) -> &mut ceo::Source;
        fn gmt(&mut self) -> &mut ceo::Gmt;
        fn detector(&mut self) -> &mut ceo::Imaging;
        fn lenslet_array(&self) -> LensletArray;
        fn pixel_scale(&mut self) -> f64;
    }

    #[derive(Copy, Clone)]
    pub struct LensletArray {
        pub n_side_lenslet: i32,
        pub lenslet_size: f64,
    }

    pub struct SH48 {
        pub gmt: ceo::Gmt,
        m1_n_mode: u64,
        pub sensor: ceo::Imaging,
        pub optics: LensletArray,
        pub n_px_lenslet: i32,
        pub n_px_framelet: i32,
        binning: i32,
        pub guide_star: ceo::Source,
        guide_star_band: String,
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
        fn through(&mut self) -> &mut Self {
            self.guide_star
                .through(&mut self.gmt)
                .xpupil()
                .through(&mut self.sensor);
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
            0.5*(self.binning as f64)*self.guide_star.wavelength()/self.optics.lenslet_size
        }
    }
    impl Drop for SH48 {
        fn drop(&mut self) {
            drop(&self.gmt);
            drop(&self.sensor);
            drop(&self.guide_star);
        }
    }

    pub struct Probe<'a, S: Sensor, T, R> {
        probe_coordinates: astrotools::SkyCoordinates,
        guide_star_magnitude: f64,
        pub sensor: &'a mut S,
        detector_integration_duration: u32,
        pub detector_frame: Vec<f32>,
        pub sensor_data0: ceo::Centroiding,
        pub sensor_data: ceo::Centroiding,
        sender: Sender<T>,
        receiver: Receiver<R>,
    }
    impl<'a, S, T, R> Probe<'a, S, T, R>
    where
        S: Sensor,
    {
        pub fn new(
            probe_coordinates: astrotools::SkyCoordinates,
            guide_star_magnitude: f64,
            sensor: &'a mut S,
            detector_integration_duration: u32,
            channel: (Sender<T>, Receiver<R>),
        ) -> Probe<S, T, R> {
            Probe {
                probe_coordinates,
                guide_star_magnitude,
                sensor,
                detector_integration_duration,
                detector_frame: vec![],
                sensor_data0: ceo::Centroiding::new(),
                sensor_data: ceo::Centroiding::new(),
                sender: channel.0,
                receiver: channel.1,
            }
        }
        pub fn init(&mut self, observation: &astrotools::Observation) {
            let (z, a) = self.probe_coordinates.local_polar(observation);
            println!("({},{})", z.to_degrees() * 60., a.to_degrees());
            self.sensor.build(
                vec![z as f32],
                vec![a as f32],
                vec![self.guide_star_magnitude as f32],
            );
            let resolution = self.sensor.detector().resolution() as usize;
            self.detector_frame = vec![0f32; resolution * resolution];
            let data_units = Some(self.sensor.pixel_scale());
            println!("pixel size: {:}",data_units.unwrap().to_degrees()*3600.);
            let data_units = Some(1f64); // TODO: fix units issue
            self.sensor_data0
                .build(self.sensor.lenslet_array().n_side_lenslet as u32, data_units);
            self.sensor_data
                .build(self.sensor.lenslet_array().n_side_lenslet as u32, data_units);
            self.sensor.guide_star().set_fwhm(3.16);
        }
        pub fn calibrate_sensor(
            &mut self,
            intensity_threshold: f64,
        ) {
            self.sensor.gmt().reset();
            self.sensor.detector().reset();
            self.sensor.through();
            self.sensor_data0.process(self.sensor.detector(), None);
            let n_valid_lenslet = self.sensor_data0.set_valid_lenslets(Some(intensity_threshold),None);
            println!(
                "# valid lenslet: {} ({:.0}%)",
                n_valid_lenslet,
                100f64
                    * n_valid_lenslet as f64
                        / (self.sensor.lenslet_array().n_side_lenslet
                        * self.sensor.lenslet_array().n_side_lenslet) as f64
            );
            self.sensor.detector().reset();
        }
        pub fn update(
            &mut self,
            observation: &astrotools::Observation,
            m1_rbm: Option<&Vec<Vec<f64>>>,
            m2_rbm: Option<&Vec<Vec<f64>>>,
        ) -> Option<Vec<f32>> {
            let (z, a) = self.probe_coordinates.local_polar(observation);
            //println!("({},{})", z.to_degrees() * 60., a.to_degrees());
            self.sensor.guide_star().update(vec![z], vec![a]);
            self.sensor.gmt().update(m1_rbm, m2_rbm);
            self.sensor.through();
            if self.sensor.detector().n_frame() == self.detector_integration_duration {
                self.sensor_data
                    .process(self.sensor.detector(), Some(&self.sensor_data0));
                self.sensor.detector().reset();
                return Some(self.sensor_data.grab().valids(Some(&self.sensor_data0.valid_lenslets)))
            }
            None
        }
    }
}
