pub mod probe {

    use crate::astrotools;
    use crate::ceo;
    use crossbeam_channel::{bounded, unbounded, Receiver, Sender};

    pub trait Sensor {
        fn build(&mut self, zenith: Vec<f32>, azimuth: Vec<f32>, magnitude: Vec<f32>);
        fn through(&mut self);
        fn guide_star(&mut self) -> &mut ceo::Source;
        fn gmt(&mut self) -> &mut ceo::Gmt;
        fn wfs(&mut self) -> &mut ceo::Imaging;
    }

    pub struct SH48 {
        gmt: ceo::Gmt,
        m1_n_mode: u64,
        sensor: ceo::Imaging,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        n_px_framelet: i32,
        binning: i32,
        lenslet_size: f64,
        pub guide_star: ceo::Source,
        guide_star_band: String,
    }
    impl SH48 {
        pub fn new() -> SH48 {
            SH48 {
                gmt: ceo::Gmt::new(),
                m1_n_mode: 17,
                sensor: ceo::Imaging::new(),
                n_side_lenslet: 48,
                n_px_lenslet: 16,
                n_px_framelet: 8,
                binning: 3,
                lenslet_size: 25.5 / 48.,
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
                self.n_side_lenslet,
                self.n_px_lenslet,
                dft_osf,
                n_px_imagelet,
                self.binning,
            );
            self.guide_star = ceo::Source::new(
                1,
                self.n_side_lenslet as f64 * self.lenslet_size,
                self.n_side_lenslet * self.n_px_lenslet + 1,
            );
            self.guide_star
                .build(&self.guide_star_band, zenith, azimuth, magnitude);
        }
        fn through(&mut self) {
            self.guide_star
                .through(&mut self.gmt)
                .xpupil()
                .through(&mut self.sensor);
        }
        fn guide_star(&mut self) -> &mut ceo::Source {
            &mut self.guide_star
        }
        fn gmt(&mut self) -> &mut ceo::Gmt {
            &mut self.gmt
        }
        fn wfs(&mut self) -> &mut ceo::Imaging {
            &mut self.sensor
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
            channel: (Sender<T>, Receiver<R>),
        ) -> Probe<S, T, R> {
            Probe {
                probe_coordinates,
                guide_star_magnitude,
                sensor,
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
        }
        pub fn update(
            &mut self,
            observation: &astrotools::Observation,
            m1_rbm: Option<Vec<Vec<f64>>>,
            m2_rbm: Option<Vec<Vec<f64>>>,
        ) {
            let (z, a) = self.probe_coordinates.local_polar(observation);
            //println!("({},{})", z.to_degrees() * 60., a.to_degrees());
            self.sensor.guide_star().update(vec![z], vec![a]);
            self.sensor.gmt().update(m1_rbm, m2_rbm);
            self.sensor.through();
            self.sensor.wfs();
        }
    }
}
