use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::agws;
use gicsdom::agws::probe::GmtState;
use gicsdom::agws::probe::Sensor;
use gicsdom::astrotools;
use gicsdom::ceo;
use gicsdom::DomeSeeing;
use ndarray::Array;
use std::f32;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

const N_PROBE: usize = 3;

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

    let telescope_channel: Vec<(Sender<GmtState>, Receiver<GmtState>)> = vec![unbounded(); N_PROBE];
    let probe_channel: Vec<(Sender<Option<Vec<f32>>>, Receiver<Option<Vec<f32>>>)> =
        vec![unbounded(); N_PROBE];
    let science_channel: (Sender<GmtState>, Receiver<GmtState>) = unbounded();

    let obs = Arc::new(Mutex::new(astrotools::Observation::from_date_utc(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        astrotools::Time::from_date_utc(2012, 6, 10, 4, 1, 3),
        astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107),
        sampling_time,
        duration,
    )));
    //    println!("OBS: {}", obs.borrow().utc.datetime());

    let mut domeseeing = DomeSeeing::new(0, 0, "cd", 12, Some(1.0 / sampling_time));
    domeseeing.list();
    let opd = Arc::new(Mutex::new(vec![0f32; 769 * 769]));

    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let mut handles = vec![];
    for k in 0..N_PROBE {
        let _probe_channel_ = probe_channel[k].clone();
        let _telescope_channel_ = telescope_channel[k].clone();
        let _obs_ = Arc::clone(&obs);
        let _opd_ = Arc::clone(&opd);
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

            // Calibration ----------------
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
            // Calibration ----------------

            let mut gopd = ceo::CuFloat::new();
            gopd.malloc(769 * 769);

            probe.sensor.detector().reset();

            loop {
                //println!("Probe #{} locking obs",k);
                {
                    let __obs = _obs_.lock().unwrap();

                    println!(
                        "Probe: {} #{:03}",
                        __obs.utc.datetime(),
                        __obs.step,
                    );

                    if __obs.ended {
                        println!("Observation has ended!");
                        break;
                    }
                    probe.update(&__obs);
                }
                //println!("Probe #{} releasing obs",k);
                //println!("Probe #{} locking opd",k);
                {
                    let mut __opd = _opd_.lock().unwrap();
                    gopd.up(&mut __opd);
                    probe.through(Some(&mut gopd));
                }
                //println!("Probe #{} releasing opd",k);
            }
        });
        handles.push(handle);
    }

    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Calibration ----------------
    let calibration: Vec<Vec<f32>> = probe_channel.iter().map(|x| x.1.recv().unwrap().unwrap()).collect();
    let n_mode = 14;
    let m = ceo::calibrations::pseudo_inverse(calibration, n_mode);
    // Calibration ----------------

    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let science_receiver = science_channel.1.clone();
    let _obs_ = Arc::clone(&obs);
    let handle = thread::spawn(move || {
        let n_src = 1usize;
        let zen = vec![0f32;n_src];
        let azi = vec![0f32;n_src];
        let pupil_size = 25.5;
        let pupil_sampling = 512;
        let m1_n_mode = 27;

        // GMT initialization
        let mut gmt = ceo::Gmt::new();
        gmt.build(m1_n_mode, None);

        let mut src = ceo::Source::new(1, pupil_size, pupil_sampling);
        src.build("V", zen, azi, vec![0.0; n_src]);

        let mut src_pssn = ceo::PSSn::new(15e-2, 25.0);
        src_pssn.build(&mut src);

        loop {
            {
                let __obs = _obs_.lock().unwrap();
                if __obs.ended {
                    println!("Observation has ended!");
                    break;
                }
            }

            let gmt_state: GmtState = science_receiver.recv().unwrap();
            gmt.update(Some(&gmt_state.m1_rbm), Some(&gmt_state.m2_rbm));
            let src_wfe_rms =src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0];
            src_pssn.peek(&mut src);
            let onaxis_pssn = src_pssn.estimates[0];
            let pssn_spatial_uniformity = src_pssn.spatial_uniformity();
            println!("SCIENCE :: WFE RMS: {:.1} nm; PSSn: {:>6.4} {:>6.4}%",
                     src_wfe_rms, onaxis_pssn, pssn_spatial_uniformity);
        }

    });
    handles.push(handle);
    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    //let stt0 = probe.sensor.guide_star().segments_gradients();
    /*
    state.m2_rbm[0][3] = 1e-6;
    state.m2_rbm[0][4] = 1e-6;
    state.m2_rbm[2][3] = 1e-6;
    state.m2_rbm[4][4] = 1e-6;
    */

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

    {
        //print!("Plant sending GMT state ...");
        telescope_channel.iter().for_each(|x| x.0.send(state.clone()).unwrap());
        science_channel.0.send(state.clone()).unwrap();
        //println!("OK");
        let mut _opd = opd.lock().unwrap();
        /*
        *_opd = domeseeing
            .next()
            .unwrap()
            .iter()
            .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
            .collect();
        */
        *_opd = vec![0f32;769*769];
    }
    obs.lock().unwrap().next();

    let now = Instant::now();
    loop {
        //print!("Plant receiving centroids ...");
        let centroids: Vec<Option<Vec<f32>>> = probe_channel.iter().map(|x| x.1.recv().unwrap()).collect();
        //println!("OK");
        if centroids.iter().all(|x| x.is_some()) {
            println!("Some data!");
            let mut cat_centroids: Vec<f32> = Vec::new();
            centroids.iter().for_each( |x| cat_centroids.append(&mut x.clone().unwrap()));
            let slopes = Array::from_shape_vec((cat_centroids.len(), 1), cat_centroids).unwrap();
            let m2_tt = m.dot(&slopes) * u * 1e6;
            println!("M2 TT:\n{:+3.2}", m2_tt.into_shape((7, 2)).unwrap());
        }

        //println!("Plant locking obs");
        {
            //print!("Plant sending GMT state ...");
            {
                telescope_channel.iter().for_each(|x| x.0.send(state.clone()).unwrap());
                science_channel.0.send(state.clone()).unwrap();
                let mut _opd = opd.lock().unwrap();
                /*
                *_opd = domeseeing
                    .next()
                    .unwrap()
                    .iter()
                    .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                    .collect();
                */
                *_opd = vec![0f32;769*769];
            }
            //println!("OK");
            let mut _obs_ = obs.lock().unwrap();
            println!("Plant: {} #{:03}", _obs_.utc.datetime(), _obs_.step,);
            if _obs_.next().is_none() {
                println!("Breaking!");
                probe_channel.iter().for_each(|x| {x.1.recv().unwrap();});
                break;
            }
        }
        //println!("Plant releasing obs");
    }
    println!("ELAPSED TIME: {}s", now.elapsed().as_secs());

    for handle in handles {
        handle.join().unwrap();
    }
}
