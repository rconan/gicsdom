use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::agws;
use gicsdom::agws::probe::GmtState;
use gicsdom::agws::probe::Sensor;
use gicsdom::astrotools;
use gicsdom::ceo;
use gicsdom::DomeSeeing;
use ndarray::Array;
use std::cell::RefCell;
use std::f32;
use std::rc::Rc;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

const N_PROBE: usize = 1;

fn main() {
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
    let mut state = GmtState::new();

    let telescope_channel: (Sender<GmtState>, Receiver<GmtState>) = bounded(0);
    let probe_channel: (Sender<Option<Vec<f32>>>, Receiver<Option<Vec<f32>>>) = bounded(0);

    let obs = Arc::new(Mutex::new(astrotools::Observation::from_date_utc(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        astrotools::Time::from_date_utc(2012, 6, 10, 4, 1, 3),
        astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107),
        sampling_time,
        duration,
    )));
    //    println!("OBS: {}", obs.borrow().utc.datetime());

    /*
    let mut domeseeing = DomeSeeing::new(0, 0, "cd", 12, Some(1.0 / sampling_time));
    domeseeing.list();
     */

    let mut handles = vec![];
    let _probe_channel_ = probe_channel.clone();
    let _telescope_channel_ = telescope_channel.clone();
    let _obs_ = Arc::clone(&obs);
    let handle = thread::spawn(move || {
        let mut sh48 = agws::probe::SH48::new();
        let optics = sh48.optics;
        let mut probe = agws::probe::Probe::new(
            astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107),
            18.,
            &mut sh48,
            sh48_integration,
            (_probe_channel_.0.clone(), _telescope_channel_.1),
        );
        {
            let __obs = _obs_.lock().unwrap();
            probe.init(&__obs);
        }
        //        probe.init(&obs.borrow());
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
        _probe_channel_.0.send(Some(calibration)).unwrap();

        let mut gopd = ceo::CuFloat::new();
        gopd.malloc(769 * 769);

        probe.sensor.detector().reset();

        loop {
            let mut opd = vec![0f32;769*769];
            /*
            let mut opd: Vec<f32> = domeseeing
                .next()
                .unwrap()
                .iter()
                .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                .collect();
            //        let opd = domeseeing.next().unwrap();
            */
            gopd.up(&mut opd);
            //println!("Probe locking obs");
            {
                let __obs = _obs_.lock().unwrap();
                /*
                println!(
                    "Probe: {} #{:03}",
                    __obs.utc.datetime(),
                    __obs.step,
                );
                */
                if __obs.ended {
                    println!("Observation has ended!");
                    break;
                }
                probe.update(&__obs);
            }
            //println!("Probe releasing obs");
            probe.through(None);//Some(&mut gopd));

        }
    });
    handles.push(handle);

    let calibration = probe_channel.1.recv().unwrap().unwrap();
    let n_mode = 14;
    let n_data = calibration.len() / n_mode;
    let m = ceo::calibrations::pseudo_inverse(calibration, n_data, n_mode);

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
    //    probe.sensor.detector().reset();

    let u: f32 = 1.5 * 0.715e-6 * 48.0 / 25.5;
    println!("U={}", 3600.0 * u.to_degrees());

    print!("Plant sending GMT state ...");
    telescope_channel.0.send(state.clone()).unwrap();
    println!("OK");

    let now = Instant::now();
    loop {

        //print!("Plant receiving centroids ...");
        let centroids = probe_channel.1.recv().unwrap();
        //println!("OK");
        if centroids.is_some() {
            println!("Some data!");
            let slopes = Array::from_shape_vec((n_data, 1), centroids.unwrap()).unwrap();
            let m2_tt = m.dot(&slopes) * u * 1e6;
            println!("M2 TT:\n{:+3.2}", m2_tt.into_shape((7, 2)).unwrap());
        }

        //println!("Plant locking obs");
        {
            //print!("Plant sending GMT state ...");
            telescope_channel.0.send(state.clone()).unwrap();
            //println!("OK");
            let mut _obs_ = obs.lock().unwrap();
            println!(
                "Plant: {} #{:03}",
                _obs_.utc.datetime(),
                _obs_.step,
            );
            if _obs_.next().is_none() {
                println!("Breaking!");
                probe_channel.1.recv().unwrap();
                break;
            }
        }
        //println!("Plant releasing obs");
    }
    println!("ELAPSED TIME: {}s", now.elapsed().as_millis());

    for handle in handles {
        handle.join().unwrap();
    }
}
