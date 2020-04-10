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
use termion;
use serde::Deserialize;
use csv::Reader;

const N_PROBE: usize = 3;

#[derive(Default, Debug, Deserialize)]
//#[serde(rename_all = "PascalCase")]
struct StarField {
    index: u32,
    datetime: String,
    ra0: f64,
    dec0: f64,
    ra1: f64,
    dec1: f64,
    v1: f64,
    p1: u8,
    ra2: f64,
    dec2: f64,
    v2: f64,
    p2: u8,
    ra3: f64,
    dec3: f64,
    v3: f64,
    p3: u8,
    ra4: f64,
    dec4: f64,
    v4: f64,
    p4: u8,
}
#[derive(Default,Debug, Deserialize)]
struct ScienceField {
    zenith: f64,
    azimuth: f64,
}


fn main() {
    println! {"{}",termion::clear::All};

    // STARFIELD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let mut field: StarField = Default::default();
    let mut rdr = Reader::from_path("guide_stars.rs.csv").unwrap();
    for result in rdr.deserialize() {
        field = result.unwrap();
//        println!("{:?}", field);
    }
    let datetime = field.datetime;
    let telescope = (field.ra0, field.dec0);
    let probes = [
        (field.ra2, field.dec2),
        (field.ra3, field.dec3),
        (field.ra4, field.dec4),
    ];
    let probe_ids = [field.p2 as i32, field.p3 as i32, field.p4 as i32];
    let probe_sh48_gs_mag = [field.v2, field.v3, field.v4];
    // STARFIELD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    println!("{}{}######## MAIN #########{}",
             termion::style::Invert,
             termion::cursor::Goto(1, 1),
             termion::style::Reset);

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
            println!(
                "{}{}###### AGWS PROBE #{} ######{}",
                termion::cursor::Goto(1, (10 * (k + 1)) as u16),
                termion::style::Invert,
                k + 1,
                termion::style::Reset
            );
            let mut sh48 = agws::probe::SH48::new();
            let optics = sh48.optics;
            let mut probe = agws::probe::Probe::new(
                astrotools::SkyCoordinates::new(probes[k].0,probes[k].1),
                probe_sh48_gs_mag[k],
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
            println!(
                "{}Calibrating ... ",
                termion::cursor::Goto(1, (10 * (k + 1) + 2) as u16)
            );
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
            println!(
                "{}Calibrated in {}s",
                termion::cursor::Goto(1, (10 * (k + 1) + 2) as u16),
                now.elapsed().as_millis()
            );
            _probe_channel_.0.send(Some(calibration)).unwrap();
            // Calibration ----------------

            let mut gopd = ceo::CuFloat::new();
            gopd.malloc(769 * 769);

            probe.sensor.detector().reset();

            loop {
                println!(
                    "{}{}OBS locked",
                    termion::cursor::Goto(1, (10 * (k + 1) + 3) as u16),
                    termion::clear::CurrentLine
                );
                {
                    let __obs = _obs_.lock().unwrap();

                    println!(
                        "{}Probe: {} #{:03}",
                        termion::cursor::Goto(1, (10 * (k + 1) + 1) as u16),
                        __obs.utc.datetime(),
                        __obs.step,
                    );

                    if __obs.ended {
                        //println!("Observation has ended!");
                        println!(
                            "{}OBS released",
                            termion::cursor::Goto(1, (10 * (k + 1) + 3) as u16)
                        );
                        break;
                    }
                    probe.update(&__obs);
                }
                //println!("Probe #{} releasing obs",k);
                println!(
                    "{}OBS released",
                    termion::cursor::Goto(1, (10 * (k + 1) + 3) as u16)
                );
                //println!("Probe #{} locking opd",k);
                println!(
                    "{}{}OPD locked",
                    termion::cursor::Goto(1, (10 * (k + 1) + 4) as u16),
                    termion::clear::CurrentLine
                );
                {
                    let mut __opd = _opd_.lock().unwrap();
                    gopd.up(&mut __opd);
                    probe.through(Some(&mut gopd));
                }
                //println!("Probe #{} releasing opd",k);
                println!(
                    "{}OPD released",
                    termion::cursor::Goto(1, (10 * (k + 1) + 4) as u16)
                );
            }
        });
        handles.push(handle);
    }

    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Calibration ----------------
    let calibration: Vec<Vec<f32>> = probe_channel
        .iter()
        .map(|x| x.1.recv().unwrap().unwrap())
        .collect();
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
        src_pssn.build(&mut src.through(&mut gmt).xpupil());

        loop {
            {
                let __obs = _obs_.lock().unwrap();
                if __obs.ended {
                    println!("{}Science Observation has ended!",
                             termion::cursor::Goto(1,7 as u16));
                    break;
                }
            }

            println!(
                "{}{}Science <<<",
                termion::cursor::Goto(1,6),
                termion::clear::CurrentLine
            );
            let gmt_state: GmtState = science_receiver.recv().unwrap();
            println!(
                "{}Science <<< OK!",
                termion::cursor::Goto(1,6)
            );

            gmt.update(Some(&gmt_state.m1_rbm), Some(&gmt_state.m2_rbm));
            let src_wfe_rms =src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0];
            src_pssn.peek(&mut src);
            let onaxis_pssn = src_pssn.estimates[0];
            let pssn_spatial_uniformity = src_pssn.spatial_uniformity();
            println!("{}SCIENCE :: WFE RMS: {:.1} nm; PSSn: {:>6.4} {:>6.4}%",
                     termion::cursor::Goto(1,5 as u16),
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
    //println!("U={}", 3600.0 * u.to_degrees());

    println!(
        "{}{}OPD locked",
        termion::cursor::Goto(1,4),
        termion::clear::CurrentLine
    );
    {
        //print!("Plant sending GMT state ...");
        telescope_channel
            .iter()
            .for_each(|x| x.0.send(state.clone()).unwrap());
        science_channel.0.send(state.clone()).unwrap();
        //println!("OK");
        let mut _opd = opd.lock().unwrap();
        *_opd = domeseeing
            .next()
            .unwrap()
            .iter()
            .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
            .collect();
//        *_opd = vec![0f32; 769 * 769];
    }
    println!(
        "{}OPD released",
        termion::cursor::Goto(1,4)
    );
    println!(
        "{}{}OBS locked",
        termion::cursor::Goto(1,3),
        termion::clear::CurrentLine
    );
    {
        obs.lock().unwrap().next();
    }
    println!(
        "{}OBS released",
        termion::cursor::Goto(1,3)
    );
 
    let now = Instant::now();
    loop {
        //print!("Plant receiving centroids ...");
        let centroids: Vec<Option<Vec<f32>>> =
            probe_channel.iter().map(|x| x.1.recv().unwrap()).collect();
        //println!("OK");
        if centroids.iter().all(|x| x.is_some()) {
            //println!("Some data!");
            let mut cat_centroids: Vec<f32> = Vec::new();
            centroids
                .iter()
                .for_each(|x| cat_centroids.append(&mut x.clone().unwrap()));
            let slopes = Array::from_shape_vec((cat_centroids.len(), 1), cat_centroids).unwrap();
            let m2_tt = m.dot(&slopes) * u * 1e6;
            println!(
                "{}M2 TT:",
                termion::cursor::Goto(1, 40)
            );
            println!(
                "{}{:+3.2}",
                termion::cursor::Goto(1, 41),
                m2_tt.into_shape((7, 2)).unwrap()
            );
        }

        //println!("Plant locking obs");
        //print!("Plant sending GMT state ...");
        println!(
            "{}{}OPD locked",
            termion::cursor::Goto(1,4),
            termion::clear::CurrentLine
        );
        {
            telescope_channel
                .iter()
                .for_each(|x| x.0.send(state.clone()).unwrap());
            science_channel.0.send(state.clone()).unwrap();
            let mut _opd = opd.lock().unwrap();
            *_opd = domeseeing
                .next()
                .unwrap()
                .iter()
                .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                .collect();
//            *_opd = vec![0f32; 769 * 769];
        }
        //println!("OK");
        println!(
            "{}OPD released",
            termion::cursor::Goto(1,4)
        );
        println!(
            "{}{}OBS locked",
            termion::cursor::Goto(1,3),
            termion::clear::CurrentLine
        );
        {
            let mut _obs_ = obs.lock().unwrap();
            println!(
                "{}Plant: {} #{:03}",
                termion::cursor::Goto(1, 2),
                _obs_.utc.datetime(),
                _obs_.step,
            );
            if _obs_.next().is_none() {
                println!(
                    "{}Plant: {} #{:03}",
                    termion::cursor::Goto(1, 2),
                    _obs_.utc.datetime(),
                    _obs_.step,
                );
                //println!("Breaking!");
                println!(
                    "{}OBS released",
                    termion::cursor::Goto(1,3)
                );
                probe_channel.iter().for_each(|x| {
                    x.1.recv().unwrap();
                });
                break;
            }
        }
        println!(
            "{}OBS released",
            termion::cursor::Goto(1,3)
        );
        //println!("Plant releasing obs");
    }
    println!(
        "{}ELAPSED TIME: {}s",
        termion::cursor::Goto(1, 39),
        now.elapsed().as_secs()
    );

    for handle in handles {
        handle.join().unwrap();
    }
    println!("{}", termion::cursor::Goto(1, 50));
}
