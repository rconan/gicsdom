use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use csv::Reader;
use gicsdom::agws;
use gicsdom::agws::probe::Sensor;
use gicsdom::agws::probe::{GmtState, Message};
use gicsdom::astrotools;
use gicsdom::ceo;
use gicsdom::ceo::calibrations::{Mirror, RigidBodyMotion};
use gicsdom::DomeSeeing;
use ndarray::Array;
use serde::Deserialize;
//use std::ops::Range;
//use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;
use std::{f32, f64};
use termion;

const N_PROBE: usize = 3;
const WITH_DOME_SEEING: bool = true;
const SIM_DURATION_MN: f64 = 10f64;

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
    let date = datetime
        .rsplit("T")
        .last()
        .unwrap()
        .split('-')
        .map(|x| x.parse::<f64>().unwrap())
        .collect::<Vec<f64>>();
    let time = datetime
        .split("T")
        .last()
        .unwrap()
        .split(':')
        .map(|x| x.parse::<f64>().unwrap())
        .collect::<Vec<f64>>();
    let telescope = (field.ra0, field.dec0);
    let probes = [
        (field.ra2, field.dec2),
        (field.ra3, field.dec3),
        (field.ra4, field.dec4),
    ];
    let _probe_ids = [field.p2 as i32, field.p3 as i32, field.p4 as i32];
    let probe_sh48_gs_mag = [field.v2, field.v3, field.v4];
    // STARFIELD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    println!(
        "{}{}TCS{}",
        termion::style::Invert,
        termion::cursor::Goto(4, 1),
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
    let duration = SIM_DURATION_MN * 60.0;
    let n_step = (duration / sampling_time) as usize;
    let sh48_integration = (30.0 / sampling_time) as u32;
    let mut state0 = GmtState::new();

    let telescope_channel: Vec<(Sender<Message>, Receiver<Message>)> = vec![bounded(0); N_PROBE];
    let probe_channel: Vec<(Sender<Option<Vec<f32>>>, Receiver<Option<Vec<f32>>>)> =
        vec![bounded(0); N_PROBE];
    let science_channel: (Sender<Message>, Receiver<Message>) = unbounded();

    let mut obs = astrotools::Observation::from_date_utc(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        astrotools::Time::from_date_utc(
            date[0] as i32,
            date[1] as u8,
            date[2] as u8,
            time[0] as u8,
            time[1] as u8,
            time[2] as u8,
        ),
        astrotools::SkyCoordinates::new(telescope.0, telescope.1),
        sampling_time,
        duration,
    );
    //    println!("OBS: {}", obs.borrow().utc.datetime());

    let mut domeseeing = DomeSeeing::new(0, 0, "cd", 12, Some(1.0 / sampling_time));
    domeseeing.list();

    let mirrors = vec![Mirror::M1,Mirror::M2];
    //    let rbms = vec![RigidBodyMotion::Rxyz(1e-6, Some(0..2))];
    let rbms = (0..7)
        .into_iter()
        .map(|k| {
            if k == 6 {
                vec![
                    RigidBodyMotion::Txyz(1e-6, Some(0..2)),
                    RigidBodyMotion::Rxyz(1e-6, Some(0..2)),
                ]
            } else {
                vec![
                    RigidBodyMotion::Txyz(1e-6, None),
                    RigidBodyMotion::Rxyz(1e-6, None),
                ]
            }
        })
        .collect::<Vec<Vec<ceo::calibrations::RigidBodyMotion>>>();
    let n_mode = mirrors.len()
        * rbms
            .iter()
            .fold(0, |a, x| a + x.iter().fold(0, |a, x| a + x.n_mode()));
    let n_threshold = Some(36);

    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let mut handles = vec![];
    for k in 0..N_PROBE {
        let _probe_channel_ = probe_channel[k].clone();
        let _telescope_channel_ = telescope_channel[k].clone();
        let _mirrors_ = mirrors.clone();
        let _rbms_ = rbms.clone();
        let handle = thread::spawn(move || {
            println!(
                "{}{}AGWS PROBE #{}{}",
                termion::cursor::Goto(4, (10 * (k + 1)) as u16),
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
                (_probe_channel_.0.clone(), None),
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
            let (z, a) = probe.coordinates.local_polar(&msg.obs);
            let calibration = m2_tt
                .build(
                    z as f32,
                    a as f32,
                    &probe.sensor_data0.valid_lenslets,
                    Some(m1_n_mode),
                    None,
                )
                .calibrate(_mirrors_, _rbms_);
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
                    termion::cursor::Goto(1, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Red),
                    termion::style::Reset
                );
                let msg = _telescope_channel_.1.recv().unwrap();
                println!(
                    "{}{}R{}",
                    termion::cursor::Goto(1, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Green),
                    termion::style::Reset
                );
                let obs = msg.obs;
                if obs.ended {
                    break;
                }
                let (z, a) = probe.update(&obs, &msg.state);
                let (alt, az, prlc) = probe.coordinates.alt_az_parallactic_deg(&obs);
                println!(
                    "{} {}  {:.1} [{:.3},{:.3}] [{:.3},{:.3},{:.3}] [{:.3},{:+.3}]",
                    termion::cursor::Goto(1, (10 * (k + 1) + 1) as u16),
                    obs.utc.datetime(),
                    probe_sh48_gs_mag[k],
                    probe.coordinates.radec.0.to_degrees(),
                    probe.coordinates.radec.1.to_degrees(),
                    alt,
                    az,
                    prlc,
                    z.to_degrees() * 60.0,
                    a.to_degrees()
                );
                println!(
                    "{}{}T{}",
                    termion::cursor::Goto(2, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Red),
                    termion::style::Reset
                );
                gopd.up(&mut msg.opd.unwrap());
                let c = probe.through(Some(&mut gopd));
                if c.is_some() {
                    println!(
                        "{}Centroids: sum={:0}",
                        termion::cursor::Goto(0, (10 * (k + 1) + 4) as u16),
                        c.clone().unwrap().into_iter().sum::<f32>()
                    );
                }
                _probe_channel_.0.send(c).unwrap();
                println!(
                    "{}{}T{}",
                    termion::cursor::Goto(2, (10 * (k + 1)) as u16),
                    termion::color::Bg(termion::color::Green),
                    termion::style::Reset
                );
            }
        });
        handles.push(handle);
    }
    // AGWS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    println!(
        "{}{}SCIENCE{}",
        termion::style::Invert,
        termion::cursor::Goto(4, 5),
        termion::style::Reset
    );
    let mut science_field: ScienceField = Default::default();
    let mut zen = vec![];
    let mut azi = vec![];
    let mut rdr = Reader::from_path("KPP_field_sampler.csv").unwrap();
    for result in rdr.deserialize() {
        science_field = result.unwrap();
        zen.push((science_field.zenith.to_radians() / 60.) as f32);
        azi.push(science_field.azimuth.to_radians() as f32);
    }

    let science_receiver = science_channel.1.clone();
    let handle = thread::spawn(move || {
        let n_src = 1 as usize; //21i32;
        zen.truncate(n_src);
        azi.truncate(n_src);
        //let zen = vec![0f32; n_src];
        //let azi = vec![0f32; n_src];
        let pupil_size = 25.5;
        let pupil_sampling = 512;
        let m1_n_mode = 27;

        // GMT initialization
        let mut gmt = ceo::Gmt::new();
        gmt.build(m1_n_mode, None);

        let mut src = ceo::Source::new(n_src as i32, pupil_size, pupil_sampling);
        src.build("V", zen, azi, vec![0.0; n_src]);

        let mut src_pssn = ceo::PSSn::new(15e-2, 25.0);
        src_pssn.build(&mut src.through(&mut gmt).xpupil());

        loop {
            let msg = science_receiver.recv().unwrap();
            let obs = msg.obs;
            if obs.ended {
                println!(
                    "{}{}Science Observation has ended!",
                    termion::cursor::Goto(1, 7),
                    termion::clear::CurrentLine,
                );
                break;
            }

            gmt.update(Some(&msg.state.m1_rbm), Some(&msg.state.m2_rbm));
            let src_wfe_rms = src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0];
            src_pssn.peek(&mut src);
            let onaxis_pssn = src_pssn.estimates[0];
            let pssn_spatial_uniformity = src_pssn.spatial_uniformity();
            println!(
                "{} {}  WFE RMS: {:>4.0}nm; PSSn: {:>6.4} {:>6.4}%",
                termion::cursor::Goto(1, 6),
                obs.utc.datetime(),
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

    let msg = Message {
        obs: obs.clone(),
        state: GmtState::new(),
        opd: None,
    };
    telescope_channel
        .iter()
        .for_each(|x| x.0.send(msg.clone()).unwrap());
    science_channel.0.send(msg.clone()).unwrap();

    // Calibration ----------------
    let calibration: Vec<Vec<f32>> = probe_channel
        .iter()
        .map(|x| x.1.recv().unwrap().unwrap())
        .collect();
    let n_data = calibration
        .iter()
        .map(|x| x.len() / n_mode)
        .collect::<Vec<usize>>();
    let m = ceo::calibrations::pseudo_inverse(calibration, n_mode, n_threshold);
    println!(
        "{} Interaction matrix: [{:?};{}]",
        termion::cursor::Goto(1, 3),
        n_data,
        n_mode
    );
    // Calibration ----------------
    /*
    state0.m2_rbm[0][3] = 1e-6;
    state0.m2_rbm[0][4] = 1e-6;
    state0.m2_rbm[2][3] = 1e-6;
    state0.m2_rbm[4][4] = 1e-6;
    state0.m2_rbm[0][3] = 1e-6;
    state0.m2_rbm[0][4] = 1e-6;
    */
    println!("{}Sending", termion::cursor::Goto(1, 2));
    let msg = Message {
        obs: obs.clone(),
        state: state0.clone(),
        opd: if WITH_DOME_SEEING {
            Some(
                domeseeing
                    .next()
                    .unwrap()
                    .iter()
                    .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                    .collect(),
            )
        } else {
            Some(vec![0f32; 769 * 769])
        },
    };
    telescope_channel
        .iter()
        .for_each(|x| x.0.send(msg.clone()).unwrap());
    science_channel.0.send(msg.clone()).unwrap();

    obs.next();

    let u: f32 = 1.5 * 0.715e-6 * 48.0 / 25.5;
    //println!("U={}", 3600.0 * u.to_degrees());
    let gain = 0.5;
    let mut state = state0.clone();
    let now = Instant::now();
    loop {
        //print!("Plant receiving centroids ...");
        println!(
            "{}{}R{}",
            termion::cursor::Goto(1, 1),
            termion::color::Bg(termion::color::Red),
            termion::style::Reset
        );
        let mut centroids: Vec<Option<Vec<f32>>> =
            probe_channel.iter().map(|x| x.1.recv().unwrap()).collect();
        if centroids.iter().all(|x| x.is_some()) {
            println!("{}", termion::cursor::Goto(1, 41));
            let mut k = 0;
            while k < N_PROBE - 1 {
                let l = centroids[k].as_ref().unwrap().len();
                if l == n_data[k] {
                    k += 1;
                } else {
                    let c = centroids.remove(k);
                    centroids.push(c);
                    k = 0;
                }
            }
            //let n_data: Vec<usize> = centroids.iter().map(|x| x.clone().unwrap().len()).collect();
            //println!("{}{:?}", termion::cursor::Goto(50, 5), n_data);
            //println!("Some data!");
            let mut cat_centroids: Vec<f32> = Vec::new();
            centroids
                .iter()
                .for_each(|x| cat_centroids.append(&mut x.clone().unwrap()));
            let slopes = Array::from_shape_vec((cat_centroids.len(), 1), cat_centroids).unwrap();
            let m2_tt = (-gain * u) * m.dot(&slopes);
            state.acc_from_array(mirrors.clone(), rbms.clone(), &m2_tt);
            println!("{}{}", termion::cursor::Goto(1, 40), state);
            /*
            println!(
                "{}{:+3.2}",
                termion::cursor::Goto(1, 41),
                m2_tt.into_shape((7, 2)).unwrap()
            );
            */
        }
        println!(
            "{}{}R{}",
            termion::cursor::Goto(1, 1),
            termion::color::Bg(termion::color::Green),
            termion::style::Reset
        );

        println!(
            "{}{}T{}",
            termion::cursor::Goto(2, 1),
            termion::color::Bg(termion::color::Red),
            termion::style::Reset
        );
        let msg = Message {
            obs: obs.clone(),
            state: state.clone(),
            opd: if WITH_DOME_SEEING {
                Some(
                    domeseeing
                        .next()
                        .unwrap()
                        .iter()
                        .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
                        .collect(),
                )
            } else {
                Some(vec![0f32; 769 * 769])
            },
        };
        telescope_channel
            .iter()
            .for_each(|x| x.0.send(msg.clone()).unwrap());
        science_channel.0.send(msg.clone()).unwrap();
        println!(
            "{}{}T{}",
            termion::cursor::Goto(2, 1),
            termion::color::Bg(termion::color::Green),
            termion::style::Reset
        );

        let (alt, az, prlc) = obs.object.alt_az_parallactic_deg(&obs);
        println!(
            "{} {}  [{:7.3},{:7.3}] [{:7.3},{:7.3},{:7.3}] #{:}/{:}",
            termion::cursor::Goto(1, 2),
            obs.utc.datetime(),
            obs.object.radec.0.to_degrees(),
            obs.object.radec.1.to_degrees(),
            alt,
            az,
            prlc,
            obs.step,
            n_step,
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
