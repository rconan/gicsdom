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
        "2012-06-10T04:01:03",
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
    );
    let t = astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107);
    let p = astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107);

    let mut hprobe = agws::probe::Handler::new(
        p,
        0.,
        agws::probe::SH48::new(),
        (probe_channel.0, telescope_channel.1),
    );
    hprobe.init(&obs, &t);
    hprobe.sensor.through();
    println!("{:?}", hprobe.sensor.get_guide_star().wfe_rms_10e(-9));
    obs.add_seconds(30.);
    println!("{}", obs.datetime);
    hprobe.update(&obs, &t, None, None);

    
    let m1_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
    let m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
    obs.add_seconds(30.);
    println!("{}", obs.datetime);
    hprobe.update(&obs, &t, Some(m1_rbm), Some(m2_rbm));

    let stt0 = hprobe.sensor.get_guide_star().segments_gradients();

    let mut obs = astrotools::Observation::new(
        "2012-06-10T04:01:03",
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
    );
    let t = astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107);
    let p = astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107);

    for k in 0..10 {
        println!("{}", obs.datetime);

        let mut m1_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
        let mut m2_rbm: Vec<Vec<f64>> = vec![vec![0.; 6]; 7];
        m1_rbm[k%7][3] = 1e-6;
        m2_rbm[k%7][4] = 1e-6;
        hprobe.update(&obs, &t, Some(m1_rbm), Some(m2_rbm));
        let mut stt = hprobe.sensor.get_guide_star().segments_gradients();

        for a in 0..2 {
            for i in 0..7 {
                stt[a][i] -= stt0[a][i];
                stt[a][i] *= 1e6;
            }
        }


        println!("{:.2?}", stt);
        obs.add_seconds(30.);
    }
}
