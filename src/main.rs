//#[macro_use]
extern crate crossbeam_channel;

use console::Term;
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use csv::Reader;
use gicsdom::ceo;
use gicsdom::ceo::GmtState;
use gicsdom::{Observation, OpticalPathToSH48, SkyCoordinates};
//use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use ndarray::{s, stack, Array, Array2, Axis};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use serde::Deserialize;
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
    /*fn main_channel(&self) -> (Sender<S>, Receiver<C>) {
        (self.server.send.clone(), self.client.recv.clone())
    }*/
}

#[derive(Debug, Deserialize)]
//#[serde(rename_all = "PascalCase")]
struct Field {
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

fn main() {
    let chat_plant_sh48: Session<GmtState, Vec<f32>> = Session::new(bounded(0), bounded(0));
    let chat_plant_science: Session<GmtState, &str> = Session::new(bounded(0), unbounded());

    let mut field = Field {
        index: 0,
        datetime: String::new(),
        ra0: 0.0,
        dec0: 0.0,
        ra1: 0.0,
        dec1: 0.0,
        v1: 0.0,
        p1: 0,
        ra2: 0.0,
        dec2: 0.0,
        v2: 0.0,
        p2: 0,
        ra3: 0.0,
        dec3: 0.0,
        v3: 0.0,
        p3: 0,
        ra4: 0.0,
        dec4: 0.0,
        v4: 0.0,
        p4: 0,
    };
    let mut rdr = Reader::from_path("guide_stars.rs.csv").unwrap();
    for result in rdr.deserialize() {
        field = result.unwrap();
        println!("{:?}", field);
    }
    let datetime = field.datetime;
    let telescope = (field.ra0, field.dec0);
    let probes = [
        (field.ra2, field.dec2),
        (field.ra3, field.dec3),
        (field.ra4, field.dec4),
    ];
    let probe_ids = (field.p2 as i32, field.p3 as i32, field.p4 as i32);
    let probe_sh48_gs_mag = (field.v2,field.v3,field.v4);

    let sampling = 30.0;

    let term = Term::buffered_stdout();

    // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let sh48_term = term.clone();
    let sh48_datetime = datetime.clone();
    let (to_plant, from_plant) = chat_plant_sh48.dual_channel();
    let h_48 = thread::spawn(move || {
        let mut sh48_0 = OpticalPathToSH48::new(&sh48_datetime, telescope, probes[0], probe_ids.0);
        sh48_0.build(probe_sh48_gs_mag.0);
        let mut sh48_1 = OpticalPathToSH48::new(&sh48_datetime, telescope, probes[1], probe_ids.1);
        sh48_1.build(probe_sh48_gs_mag.1);
        let mut sh48_2 = OpticalPathToSH48::new(&sh48_datetime, telescope, probes[2], probe_ids.2);
        sh48_2.build(probe_sh48_gs_mag.2);

        loop {
            let gstate = match from_plant.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            let (z0, a0, p0) = sh48_0.update(sampling, &gstate).local();
            sh48_0.sensor.process();

            let (z1, a1, p1) = sh48_1.update(sampling, &gstate).local();
            sh48_1.sensor.process();

            let (z2, a2, p2) = sh48_2.update(sampling, &gstate).local();
            sh48_2.sensor.process();

            sh48_term
                .write_line(&format!(
                    "SH48 #{}  {:>4.2}' {:>+7.2} {:>+6.2}\nSH48 #{}  {:>4.2}' {:>+7.2} {:>+6.2}\nSH48 #{}  {:>4.2}' {:>+7.2} {:>+6.2}",
                    sh48_0.probe_id, z0, a0,p0,
                        sh48_1.probe_id, z1, a1,p1,
                        sh48_2.probe_id, z2, a2,p2
                ))
                .unwrap();

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
    let sci_term = term.clone();
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
            "WFE RMS: {:.1} nm; PSSn: {:>6.4}",
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

            sci_term
                .write_line(&format!(
                    "WFE RMS: {:.1} nm; PSSn: {:>6.4}",
                    src_wfe_rms, onaxis_pssn
                ))
                .unwrap();
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
    /*gstate0.rbm[(0, 0)] = 1e-5;
    gstate0.rbm[(8, 4)] = 1e-6;
    gstate0.rbm[(2, 3)] = -1e-6;
    gstate0.rbm[(11, 1)] = -1e-5;
    gstate0.bm[[0, 0]] = 1e-5;
    gstate0.bm[[1, 0]] = 1e-5;
    gstate0.bm[[2, 1]] = 1e-5;*/

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
    let mut sh48_0 = OpticalPathToSH48::new(&datetime, telescope, probes[0], 0);
    let _d_0 = sh48_0.build(probe_sh48_gs_mag.0).calibrate(None);
    let mut sh48_1 = OpticalPathToSH48::new(&datetime, telescope, probes[1], 1);
    let _d_1 = sh48_1.build(probe_sh48_gs_mag.1).calibrate(None);
    let mut sh48_2 = OpticalPathToSH48::new(&datetime, telescope, probes[2], 2);
    let _d_2 = sh48_2.build(probe_sh48_gs_mag.2).calibrate(None);
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

    let mut obs = Observation::new(&datetime, -29.04, -70.682);
    let pointing = SkyCoordinates::new(telescope.0, telescope.1);

    let from_sh48 = chat_plant_sh48.client.recv; // sh48_plan_chat.1;
    term.move_cursor_down(3).unwrap();
    term.write_line(&format!("------------------\n",)).unwrap();
    for _k in 0..30 {
        let (alt, az, pa) = pointing.alt_az_parallactic_deg(&obs);
        term.write_line(&format!(
            "{} {:+>6.2} {:+>6.2} {:+>6.2}",
            obs.datetime, alt, az, pa
        ))
        .unwrap();
        term.flush().unwrap();

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

        obs.add_seconds(sampling);
        //thread::sleep(Duration::from_secs(2));
        term.move_cursor_up(5).unwrap();
    }

    term.move_cursor_down(5).unwrap();
    term.write_line(&format!("------------------\n",)).unwrap();
    term.flush().unwrap();

    drop(to_science);
    drop(to_sh48);

    h_48.join().unwrap();
    h_science.join().unwrap();
}
