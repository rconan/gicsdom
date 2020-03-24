use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::agws;
use gicsdom::agws::probe::Sensor;
use gicsdom::astrotools;
use std::boxed::Box;

fn main() {
    let mut sh48 = Box::new(agws::probe::SH48::new());
    sh48.build(vec![0.], vec![0.], vec![0.]);
    sh48.through();
    println!("{:?}", sh48.guide_star.wfe_rms_10e(-9));

    let telescope_channel: (Sender<f64>, Receiver<f64>) = bounded(0);
    let probe_channel: (Sender<u64>, Receiver<u64>) = bounded(0);

    /*
    An AGWS is defined by a date and time of observation and
    the RA and DEC sky coordinates of:
     * the telescope pointing target
     * the probes guide stars
     * the guide star magnitude
     */

    let mut obs = astrotools::Observation::new(
        "2020-03-20T00:00:00",
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
    );
    let t = astrotools::SkyCoordinates::new(0., 0.);
    let p = astrotools::SkyCoordinates::new(0., 0.);

    let mut hprobe = agws::probe::Handler::new(
        obs,
        t,
        p,
        0.,
        agws::probe::SH48::new(),
        (probe_channel.0, telescope_channel.1),
    );
    hprobe.init();
    hprobe.sensor.through();
    println!("{:?}", hprobe.sensor.guide_star.wfe_rms_10e(-9));
}
