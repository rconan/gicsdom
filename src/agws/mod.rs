pub mod probe {

    use crate::astrotools;
    use crate::ceo;
    use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
    use std::any::Any;
    use std::boxed::Box;
    use std::rc::Rc;

    pub trait Sensor {
        fn build(&mut self, zenith: Vec<f32>, azimuth: Vec<f32>, magnitude: Vec<f32>);
        fn through(&mut self);
        fn get_guide_star(&mut self) ->&mut ceo::Source;
        fn get_gmt(&mut self) -> &mut ceo::Gmt;
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
        fn get_guide_star(&mut self) -> &mut ceo::Source {
            &mut self.guide_star
        }
        fn get_gmt(&mut self) -> &mut ceo::Gmt {
            &mut self.gmt
        }
    }
    impl Drop for SH48 {
        fn drop(&mut self) {
            drop(&self.gmt);
            drop(&self.sensor);
            drop(&self.guide_star);
        }
    }

    pub struct Handler<S, R> {
        probe_coordinates: astrotools::SkyCoordinates,
        guide_star_magnitude: f64,
        pub sensor: SH48,
        sender: Sender<S>,
        receiver: Receiver<R>,
    }
    impl<S, R> Handler<S, R> {
        pub fn new(
            probe_coordinates: astrotools::SkyCoordinates,
            guide_star_magnitude: f64,
            sensor: SH48,
            channel: (Sender<S>, Receiver<R>),
        ) -> Self {
            Handler {
                probe_coordinates,
                guide_star_magnitude,
                sensor,
                sender: channel.0,
                receiver: channel.1,
            }
        }
        /*
        pub fn init(
            &mut self,
            observation: &astrotools::Observation,
            telescope_coordinates: &astrotools::SkyCoordinates,
        ) {
            let (z, a) = self
                .probe_coordinates
                .local_polar(telescope_coordinates, observation);
            println!("({},{})",z.to_degrees()*60.,a.to_degrees());
            self.sensor.build(
                vec![z as f32],
                vec![a as f32],
                vec![self.guide_star_magnitude as f32],
            );
        }
        pub fn update(
            &mut self,
            observation: &astrotools::Observation,
            telescope_coordinates: &astrotools::SkyCoordinates,
            m1_rbm: Option<Vec<Vec<f64>>>,m2_rbm: Option<Vec<Vec<f64>>>
        ) {
            let (z, a) = self
                .probe_coordinates
                .local_polar(telescope_coordinates, observation);
            println!("({},{})",z.to_degrees()*60.,a.to_degrees());
            self.sensor.get_guide_star().update(vec![z], vec![a]);
            if m1_rbm.is_some() {
                for (k,rbm) in m1_rbm.unwrap().iter().enumerate() {
                    let gmt = self.sensor.get_gmt();
                    gmt.set_m1_segment_state((k+1)as i32, &rbm[..3], &rbm[3..]);
                }
            }
            if m2_rbm.is_some() {
                for (k,rbm) in m2_rbm.unwrap().iter().enumerate() {
                    let gmt = self.sensor.get_gmt();
                    gmt.set_m2_segment_state((k+1)as i32, &rbm[..3], &rbm[3..]);
                }
            }
            self.sensor.through();
        }
        */
    }
}
