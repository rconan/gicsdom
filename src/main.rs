//#[macro_use]
extern crate crossbeam_channel;

use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::ceo;
use gicsdom::ceo::GmtState;
use gicsdom::{Observation, OpticalPathToSH48, Rotation, SkyCoordinates};
use hifitime::Epoch;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use ndarray::{s, stack, Array, Array2, Axis};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::thread;
use std::time::{Duration, Instant};
use std::{f32, f64};

struct Dyad<T> {
    send: Sender<T>,
    recv: Receiver<T>,
}
struct Session<S, C> {
    server: Dyad<S>,
    client: Dyad<C>,
}
impl<S, C> Session<S, C> {
    fn new(server_: (Sender<S>, Receiver<S>), client_: (Sender<C>, Receiver<C>)) -> Self {
        Session {
            server: Dyad {
                send: server_.0,
                recv: server_.1,
            },
            client: Dyad {
                send: client_.0,
                recv: client_.1,
            },
        }
    }
    fn dual_channel(&self) -> (Sender<C>, Receiver<S>) {
        (self.client.send.clone(), self.server.recv.clone())
    }
    fn main_channel(&self) -> (Sender<S>, Receiver<C>) {
        (self.server.send.clone(), self.client.recv.clone())
    }
}

fn main() {
    let chat_plant_sh48: Session<GmtState, Vec<f32>> = Session::new(bounded(0), bounded(0));
    let chat_plant_science: Session<GmtState, &str> = Session::new(bounded(0), unbounded());

    let observation = Epoch::from_gregorian_utc_str("1998-08-10T23:10:00").unwrap();
    let jde: f64 = observation.as_jde_tt_days();
    let jday: f64 = jde - 2451545.0;
    println!("JULIAN DAY: {}", jday);
    let (_y, _m, _d, h, min, sec) = observation.as_gregorian_utc();
    let long = -1.9166667;
    let ut: f64 = (h as f64) + ((min as f64) + (sec as f64) / 60.0) / 60.0;
    let lst = 100.46 + 0.985647 * jday + long + 15.0 * ut;
    println!("LST: {}", lst);
    let mut obs = Observation::new("1998-08-10T23:10:00", 52.5, -1.9166667);
    let target = SkyCoordinates::new(16.695 * 15.0, 36.466667);
    println!("LST: {}", obs.local_sidereal_time_deg());
    println!("AltAz: {:?}", target.altaz_deg(&obs));

    println!("Date: {}", obs.datetime);
    println!("Date: {}", obs.add_seconds(30.0).datetime);
    println!("AltAz: {:?}", target.altaz_deg(&obs));

    let rot_x = Rotation::new(f64::consts::FRAC_PI_4, 2);
    println!("R: {:?}", rot_x.mat);
    let v = vec![1f64; 3];
    println!("Rotated vector: {:?}", rot_x.apply(v));

    let z = 6. * f32::consts::PI / 180. / 60.;
    let a = 2. * f32::consts::PI / 3.;

    // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let (to_plant, from_plant) = chat_plant_sh48.dual_channel();
    let h_48 = thread::spawn(move || {
        let mut sh48_0 = OpticalPathToSH48::new(1);
        sh48_0.build(vec![z], vec![0.0 * a]);
        let mut sh48_1 = OpticalPathToSH48::new(1);
        sh48_1.build(vec![z], vec![1.0 * a]);
        let mut sh48_2 = OpticalPathToSH48::new(1);
        sh48_2.build(vec![z], vec![2.0 * a]);

        loop {
            let gstate = match from_plant.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            sh48_0.gmt.update(&gstate);
            sh48_0.propagate_src();
            sh48_0.sensor.process();

            sh48_1.gmt.update(&gstate);
            sh48_1.propagate_src();
            sh48_1.sensor.process();

            sh48_2.gmt.update(&gstate);
            sh48_2.propagate_src();
            sh48_2.sensor.process();

            let mut q: Vec<f32> = Vec::with_capacity(
                (sh48_0.sensor.n_centroids + sh48_1.sensor.n_centroids + sh48_2.sensor.n_centroids)
                    as usize,
            );
            q.append(&mut sh48_0.sensor.centroids.clone());
            q.append(&mut sh48_1.sensor.centroids.clone());
            q.append(&mut sh48_2.sensor.centroids.clone());
            to_plant.send(q).unwrap();
        }
    });
    // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let (_to_plant, from_plant) = chat_plant_science.dual_channel();
    let h_science = thread::spawn(move || {
        let pupil_size = 25.5;
        let pupil_sampling = 512;
        let m1_n_mode = 27;

        // GMT initialization
        let mut gmt = ceo::Gmt::new(m1_n_mode, None);
        gmt.build();

        // Science initialization
        let mut src = ceo::Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0.0], vec![0.0], vec![0.0]);

        let mut src_wfe_rms = src.through(&mut gmt).wfe_rms_10e(-9)[0];
        let mut onaxis_pssn: f32;
        let mut src_pssn = ceo::PSSn::new(15e-2, 25.0, 0.0);
        src_pssn.build(&mut src);
        println!(
            "WFE RMS: {:.3} nm; PSSn: {}",
            src_wfe_rms,
            src_pssn.reset(&mut src)
        );

        loop {
            let gstate = match from_plant.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            gmt.update(&gstate);

            src_wfe_rms = src.through(&mut gmt).wfe_rms_10e(-9)[0];
            onaxis_pssn = src_pssn.reset(&mut src);

            println!("WFE RMS: {:.3} nm; PSSn: {}", src_wfe_rms, onaxis_pssn);
        }
    });
    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // PLANT

    // GMT initial conditions
    let m1_n_mode = 27;
    let mut gstate0 = GmtState {
        rbm: Array2::<f32>::zeros((14, 6)),
        bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
    };
    gstate0.rbm[(0, 0)] = 1e-5;
    gstate0.rbm[(8, 4)] = 1e-6;
    gstate0.rbm[(2, 3)] = -1e-6;
    gstate0.rbm[(11, 1)] = -1e-5;
    gstate0.bm[[0, 0]] = 1e-5;
    gstate0.bm[[1, 0]] = 1e-5;
    gstate0.bm[[2, 1]] = 1e-5;

    let n_rbm = 84 as usize;
    let m1_n_mode = 27;
    let n_side_lenslet = 48i32;
    let n_sensor = 3;
    let n_c: usize = (n_side_lenslet.pow(2) as usize) * 2 * n_sensor as usize;
    let gain = 0.5_f32;
    let mut _c_e: Array2<f32>;
    let mut rbm = Array2::<f32>::zeros((14, 6));
    let mut bm = Array2::<f32>::zeros((7, m1_n_mode as usize));

    let to_science = chat_plant_science.server.send; // plant_to_science.0;
    let to_sh48 = chat_plant_sh48.server.send; // plant_to_sh48.0;
    let mut gstate = GmtState {
        rbm: Array2::<f32>::zeros((14, 6)),
        bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
    };
    gstate.rbm = gstate0.rbm.clone();
    gstate.bm = gstate0.bm.clone();

    // ---------------------------------------------------
    // SYSTEM CALIBRATION
    //let sh48_from_calibrate = bounded(0);
    //let sp = sh48_from_calibrate.0;
    //thread::spawn(move || {
    //let mut sh48 = OpticalPathToSH48::new(3);
    //let mut _d = sh48.build(vec![z, z, z], vec![0.0 * a, a, 2.0 * a]).calibrate(None);
    print!("SH48 calibration");
    let now = Instant::now();
    let mut sh48_0 = OpticalPathToSH48::new(1);
    let _d_0 = sh48_0.build(vec![z], vec![0.0 * a]).calibrate(None);
    let mut sh48_1 = OpticalPathToSH48::new(1);
    let _d_1 = sh48_1.build(vec![z], vec![a]).calibrate(None);
    let mut sh48_2 = OpticalPathToSH48::new(1);
    let _d_2 = sh48_2.build(vec![z], vec![2.0 * a]).calibrate(None);
    let mut _d = stack(Axis(0), &[_d_0.view(), _d_1.view(), _d_2.view()]).unwrap();
    println!(" in {}ms", now.elapsed().as_millis());

    print!("SVD decomposition");
    let now = Instant::now();
    let n_sv = n_rbm + m1_n_mode as usize * 7;
    let (u, sig, v_t) = _d.svddc_inplace(UVTFlag::Some).unwrap();
    println!(" in {}ms", now.elapsed().as_millis());
    //    println!("Singular values:\n{}", sig);
    let mut i_sig = sig.mapv(|x| 1.0 / x);
    for k in 0..14 {
        i_sig[n_sv - k - 1] = 0.0;
    }

    let _u = u.unwrap();
    let _vt = v_t.unwrap();

    let l_sv = Array2::from_diag(&i_sig);
    print!("Computing the pseudo-inverse");
    let now = Instant::now();
    let __m: Array2<f32> = _vt.t().dot(&l_sv.dot(&_u.t()));
    println!(" in {}ms", now.elapsed().as_millis());

    //    sp.send(__m).unwrap();
    //});
    //let rp = sh48_from_calibrate.1;
    //println!("{:?}", rp.recv());
    //let __m = rp.recv().unwrap();
    // ---------------------------------------------------

    let mut obs = Observation::new("1998-08-10T23:10:00", -29.049, -70.682);
    let sampling = 30.0;

    let from_sh48 = chat_plant_sh48.client.recv; // sh48_plan_chat.1;
    for k in 0..20 {
        obs.add_seconds(k as f64 * sampling);
        println!("{}", obs.datetime);

        to_science.send(gstate.clone()).unwrap();
        to_sh48.send(gstate.clone()).unwrap();

        let centroids = match from_sh48.recv() {
            Ok(state) => state,
            Err(_) => break,
        };
        let slopes = Array::from_shape_vec((n_c, 1), centroids).unwrap();
        _c_e = __m.dot(&slopes);
        let q = gain
            * _c_e
                .slice(s![..84, ..])
                .to_owned()
                .into_shape((14, 6))
                .unwrap();
        rbm -= &q.view();
        if m1_n_mode > 0 {
            let q = gain
                * _c_e
                    .slice(s![84.., ..])
                    .to_owned()
                    .into_shape((7, m1_n_mode as usize))
                    .unwrap();
            bm -= &q.view();
        }

        gstate.rbm = gstate0.rbm.clone() + &rbm;
        gstate.bm = gstate0.bm.clone() + &bm;
    }

    drop(to_science);
    drop(to_sh48);

    h_48.join().unwrap();
    h_science.join().unwrap();
}
