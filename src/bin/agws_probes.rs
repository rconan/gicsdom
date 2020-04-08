use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::agws;
use gicsdom::agws::probe::GmtState;
use gicsdom::agws::probe::Sensor;
use gicsdom::astrotools;
use gicsdom::ceo;
use gicsdom::DomeSeeing;
use ndarray::{Array, Array2, ShapeBuilder};
use std::cell::RefCell;
use std::cmp::Ord;
use std::rc::Rc;
use std::time::Instant;
use std::{f32, f64};

fn main() {
    let telescope_channel: (Sender<GmtState>, Receiver<GmtState>) = bounded(1);
    let probe_channel: (Sender<Option<Vec<f32>>>, Receiver<Option<Vec<f32>>>) = bounded(1);

    /*
    An AGWS is defined by a date and time of observation and
    the RA and DEC sky coordinates of:
     * the telescope pointing target
     * the probes guide stars
     * the guide star magnitude
     */

    let sampling_time = 1.0;
    let duration = 30.0;
    let sh48_integration = (30.0 / sampling_time) as u32;

    let mut domeseeing = DomeSeeing::new(0, 0, "cd", 12, Some(1.0 / sampling_time));
    domeseeing.list();

    let obs = RefCell::new(astrotools::Observation::from_date_utc(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        astrotools::Time::from_date_utc(2012, 6, 10, 4, 1, 3),
        astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107),
        sampling_time,
        duration,
    ));
    println!("OBS: {}", obs.borrow().utc.datetime());

    let mut sh48 = agws::probe::SH48::new();
    let optics = sh48.optics;
    let mut probe = agws::probe::Probe::new(
        astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107),
        18.,
        &mut sh48,
        sh48_integration,
        (probe_channel.0, telescope_channel.1),
    );
    probe.init(&obs.borrow());
    probe.calibrate_sensor(0.9);

    print!("Calibrating ");
    let m1_n_mode = 17;
    let now = Instant::now();
    let mut m2_tt = ceo::Calibration::new(optics, None);
    let calibration = m2_tt
        .build(&probe.sensor_data0.valid_lenslets, Some(m1_n_mode), None)
        .calibrate(
            vec![ceo::calibrations::Mirror::M2],
            vec![
                ceo::calibrations::RigidBodyMotion::Rxyz(0, 1e-6),
                ceo::calibrations::RigidBodyMotion::Rxyz(1, 1e-6),
            ],
        );
    println!(" in {}s", now.elapsed().as_millis());
    let m = ceo::calibrations::pseudo_inverse(calibration, m2_tt.n_data, m2_tt.n_mode);

    let mut state = GmtState::new();

    //let stt0 = probe.sensor.guide_star().segments_gradients();
    state.m2_rbm[0][3] = 1e-6;
    state.m2_rbm[0][4] = 1e-6;
    state.m2_rbm[2][3] = 1e-6;
    state.m2_rbm[4][4] = 1e-6;

    //probe.update(&obs.borrow(), Some(&state.m1_rbm), Some(&state.m2_rbm));

    /*
    let mut stt = probe.sensor.guide_star().segments_gradients();
    for a in 0..2 {
        for i in 0..7 {
            stt[a][i] -= stt0[a][i];
            stt[a][i] *= 1e6;
        }
    }
    println!("{:+.2?}", stt);
*/
    probe.sensor.detector().reset();

    let u: f32 = 1.5 * 0.715e-6 * 48.0 / 25.5;
    println!("U={}", 3600.0 * u.to_degrees());

    let mut gopd = ceo::CuFloat::new();
    gopd.malloc(769 * 769);

    let now = Instant::now();
    loop {
        println!(
            "{} #{:03}: frame count: {:03}",
            obs.borrow().utc.datetime(),
            obs.borrow().step,
            probe.sensor.detector().n_frame(),
        );

        /*
        let opd = domeseeing.next().unwrap();
        let mut opd_f32: Vec<f32> = opd
            .iter()
            .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
            .collect();
        gopd.up(&mut opd_f32);
         */
        
        telescope_channel.0.send(state.clone()).unwrap();
        probe
            .update(&obs.borrow())
            .through(None);//Some(&mut gopd));

        let centroids = probe_channel.1.recv().unwrap();
        if centroids.is_some() {
            println!("Some data!");
            let slopes =
                Array::from_shape_vec((m2_tt.n_data as usize, 1), centroids.unwrap()).unwrap();
            let m2_tt = m.dot(&slopes) * u * 1e6;
            println!("M2 TT:\n{:+3.2}", m2_tt.into_shape((7, 2)).unwrap());
        }

        if obs.borrow_mut().next().is_none() {
            break;
        }
    }
    println!("ELAPSED TIME: {}s", now.elapsed().as_secs());
}
