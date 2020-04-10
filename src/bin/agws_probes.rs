use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use csv::Reader;
use gicsdom::agws;
use gicsdom::agws::probe::{GmtState,Message};
use gicsdom::agws::probe::Sensor;
use gicsdom::astrotools;
use gicsdom::ceo;
use gicsdom::DomeSeeing;
use ndarray::Array;
use serde::Deserialize;
use std::f32;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;
use termion;

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
#[derive(Default, Debug, Deserialize)]
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

    println!(
        "{}{}######## MAIN #########{}",
        termion::style::Invert,
        termion::cursor::Goto(1, 1),
        termion::style::Reset
    );

    /*
    An AGWS is defined by a date and time of observation and
    the RA and DEC sky coordinates of:
     * the telescope pointing target
     * the probes guide stars
     * the guide star magnitude
     */

    let sampling_time = 1.0;
    let duration = 30.0;
    let n_step = (duration/sampling_time) as usize;
    let sh48_integration = (30.0 / sampling_time) as u32;
    let mut state0 = GmtState::new();

    let telescope_channel: Vec<(Sender<Message>, Receiver<Message>)> = vec![bounded(0); N_PROBE];
    let probe_channel: Vec<(Sender<Option<Vec<f32>>>, Receiver<Option<Vec<f32>>>)> =
        vec![bounded(1); N_PROBE];
    let science_channel: (Sender<Message>, Receiver<Message>) = unbounded();

    let mut obs = astrotools::Observation::from_date_utc(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        astrotools::Time::from_date_utc(2012, 6, 10, 4, 1, 3),
        astrotools::SkyCoordinates::new(266.62156258190726, -27.776114821065107),
        sampling_time,
        duration,
    );
    //    println!("OBS: {}", obs.borrow().utc.datetime());

    let mut domeseeing = DomeSeeing::new(0, 0, "cd", 12, Some(1.0 / sampling_time));
    domeseeing.list();
    let opd = Arc::new(Mutex::new(vec![0f32; 769 * 769]));

    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let mut handles = vec![];
    for k in 0..N_PROBE {
        let _probe_channel_ = probe_channel[k].clone();
        let _telescope_channel_ = telescope_channel[k].clone();
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
                astrotools::SkyCoordinates::new(probes[k].0, probes[k].1),
                probe_sh48_gs_mag[k],
                &mut sh48,
                sh48_integration,
                (_probe_channel_.0.clone(),None),
            );

            let msg = _telescope_channel_.1.recv().unwrap();
            probe.init(&msg.obs);
            //        probe.init(&obs.borrow());
            let n_valid_lenslet = probe.calibrate_sensor(0.9);
            println!(
                "{}Sensor calibrated: valid lenset: {} ({:.0}%)",
                termion::cursor::Goto(1, (10 * (k + 1) + 2) as u16),
                n_valid_lenslet,
                100f64 * n_valid_lenslet as f64
                    / (optics.n_side_lenslet * optics.n_side_lenslet) as f64
            );

            // Calibration ----------------
            println!(
                "{}Interaction matrix ... ",
                termion::cursor::Goto(1, (10 * (k + 1) + 3) as u16)
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
                "{}Interaction matrix computed in {}ms",
                termion::cursor::Goto(1, (10 * (k + 1) + 3) as u16),
                now.elapsed().as_millis()
            );
            _probe_channel_.0.send(Some(calibration)).unwrap();
            // Calibration ----------------

            let mut gopd = ceo::CuFloat::new();
            gopd.malloc(769 * 769);

            probe.sensor.detector().reset();

            loop {
                println!(
                    "{}{}R{}",
                    termion::cursor::Goto(35, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Red),
                    termion::style::Reset
                );
                let msg = _telescope_channel_.1.recv().unwrap();
                println!(
                    "{}{}R{}",
                    termion::cursor::Goto(35, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Green),
                    termion::style::Reset
                );
                let obs = msg.obs;
                if obs.ended {
                    break;
                }
                let (z, a) = probe.update(&obs,&msg.state);
                println!(
                    "{} {} [{:7.3},{:7.3}]  [{:7.3},{:+8.3}]",
                    termion::cursor::Goto(1, (10 * (k + 1) + 1) as u16),
                    obs.utc.datetime(),
                    probe.coordinates.radec.0.to_degrees(),
                    probe.coordinates.radec.1.to_degrees(),
                    z.to_degrees() * 60.0,
                    a.to_degrees()
                );
                println!(
                    "{}{}T{}",
                    termion::cursor::Goto(40, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Red),
                    termion::style::Reset
                );
                gopd.up(&mut msg.opd.unwrap());
                probe.through(Some(&mut gopd));
                println!(
                    "{}{}T{}",
                    termion::cursor::Goto(40, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Green),
                    termion::style::Reset
                );
            }
        });
        handles.push(handle);
    }

    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let science_receiver = science_channel.1.clone();
    let handle = thread::spawn(move || {
        let n_src = 1usize;
        let zen = vec![0f32; n_src];
        let azi = vec![0f32; n_src];
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
            let msg =  science_receiver.recv().unwrap();
            let obs = msg.obs;
            if obs.ended {
                println!(
                    "{}Science Observation has ended!",
                    termion::cursor::Goto(1, 4 as u16)
                );
                break;
            }

            gmt.update(Some(&msg.state.m1_rbm), Some(&msg.state.m2_rbm));
            let src_wfe_rms = src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0];
            src_pssn.peek(&mut src);
            let onaxis_pssn = src_pssn.estimates[0];
            let pssn_spatial_uniformity = src_pssn.spatial_uniformity();
            println!(
                "{}SCIENCE :: WFE RMS: {:.1} nm; PSSn: {:>6.4} {:>6.4}%",
                termion::cursor::Goto(1, 3),
                src_wfe_rms,
                onaxis_pssn,
                pssn_spatial_uniformity
            );
        }
    });
    handles.push(handle);
    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    //let stt0 = probe.sensor.guide_star().segments_gradients();

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
        "{}Sending",
        termion::cursor::Goto(1, 2));
    let msg = Message {
        obs: obs.clone(),
        state: GmtState::new(),
        opd: None,
    };
    telescope_channel
        .iter()
        .for_each(|x| x.0.send(msg.clone()).unwrap());
    science_channel.0.send(msg.clone()).unwrap();
    println!(
        "{}Sending OK",
        termion::cursor::Goto(1, 2));

    // Calibration ----------------
    let calibration: Vec<Vec<f32>> = probe_channel
        .iter()
        .map(|x| x.1.recv().unwrap().unwrap())
        .collect();
    let n_mode = 14;
    let n_data = calibration.iter().map(|x|x.len()/n_mode).collect::<Vec<usize>>();
    let m = ceo::calibrations::pseudo_inverse(calibration, n_mode);
    println!(
        "{} Calibration data: {:?}",
        termion::cursor::Goto(1, 5),
        n_data);
    // Calibration ----------------

    state0.m2_rbm[0][3] = 1e-6;
    state0.m2_rbm[0][4] = 1e-6;
    state0.m2_rbm[2][3] = 1e-6;
    state0.m2_rbm[4][4] = 1e-6;
    println!(
        "{}Sending",
        termion::cursor::Goto(1, 2));
    let msg = Message {
        obs: obs.clone(),
        state: state0.clone(),
        opd: Some(domeseeing
                  .next()
                  .unwrap()
                  .iter()
                  .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                  .collect())
//        opd: Some(vec![0f32; 769 * 769]),
    };
    telescope_channel
        .iter()
        .for_each(|x| x.0.send(msg.clone()).unwrap());
    science_channel.0.send(msg.clone()).unwrap();
    println!(
        "{}Sending OK",
        termion::cursor::Goto(1, 2));
    obs.next();

    let gain = 0.5;
    let state = state0.clone();
    let now = Instant::now();
    loop {
        //print!("Plant receiving centroids ...");
        println!(
            "{}{}R{}",
            termion::cursor::Goto(35, 1),
            termion::color::Bg(termion::color::Red),
            termion::style::Reset
        );
        let mut centroids: Vec<Option<Vec<f32>>> =
            probe_channel.iter().map(|x| x.1.recv().unwrap()).collect();
        println!(
            "{}{}R{}",
            termion::cursor::Goto(35, 1),
            termion::color::Bg(termion::color::Green),
            termion::style::Reset
        );

        if centroids.iter().all(|x| x.is_some()) {
            let mut k = 0;
            while k<N_PROBE-1 {
                let l = centroids[k].as_ref().unwrap().len();
                if l==n_data[k] {
                    k += 1;
                } else {
                    let c = centroids.remove(k);
                    centroids.push( c );
                    k = 0;
                }
            }
            let n_data: Vec<usize> = centroids.iter().map(|x| x.clone().unwrap().len()).collect();
            println!(
                "{}{:?}",
                termion::cursor::Goto(50, 5),
                n_data);
            //println!("Some data!");
            let mut cat_centroids: Vec<f32> = Vec::new();
            centroids
                .iter()
                .for_each(|x| cat_centroids.append(&mut x.clone().unwrap()));
            let slopes = Array::from_shape_vec((cat_centroids.len(), 1), cat_centroids).unwrap();
            let m2_tt = m.dot(&slopes) * u * 1e6;
            println!("{}M2 TT:", termion::cursor::Goto(1, 40));
            println!(
                "{}{:+3.2}",
                termion::cursor::Goto(1, 41),
                m2_tt.into_shape((7, 2)).unwrap()
            );
        }

        println!(
            "{}{}T{}",
            termion::cursor::Goto(40, 1),
            termion::color::Bg(termion::color::Red),
            termion::style::Reset
        );
        let msg = Message {
            obs: obs.clone(),
            state: state0.clone(),
            opd: Some(domeseeing
                      .next()
                      .unwrap()
                      .iter()
                      .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                      .collect()),
            //            opd: Some(vec![0f32; 769 * 769]),
        };
        telescope_channel
            .iter()
            .for_each(|x| x.0.send(msg.clone()).unwrap());
        science_channel.0.send(msg.clone()).unwrap();
        println!(
            "{}{}T{}",
            termion::cursor::Goto(40, 1),
            termion::color::Bg(termion::color::Green),
            termion::style::Reset
        );

        println!(
                "{} {} #{:06}/{:06}",
                termion::cursor::Goto(1, 2),
                obs.utc.datetime(),
                obs.step, n_step,
        );
        if obs.ended {
            break;
        }
        obs.next();
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
