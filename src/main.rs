//#[macro_use]
extern crate crossbeam_channel;

use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::ceo;
use gicsdom::ceo::GmtState;
//use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{s,Array, Array2};
use std::f32;
use std::thread;
use std::time::{Duration};//, Instant};

use gicsdom::OpticalPathToSH48;


fn main() {
    // Create a zero-capacity channel.
    let plant_to_sh48: (Sender<GmtState>, Receiver<GmtState>) = bounded(0);
    let sh48_plan_chat: (Sender<GmtState>, Receiver<GmtState>) = bounded(0);
    let plant_to_science: (Sender<GmtState>, Receiver<GmtState>) = unbounded();
    let science_plan_chat = bounded(0);

    // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let r = plant_to_sh48.1;
    let s = sh48_plan_chat.0;
    thread::spawn(move || {

        let mut sh48 = OpticalPathToSH48::new();
        sh48.build();

        // ---------------------------------------------------
        // SYSTEM CALIBRATION
        let sh48_from_calibrate = bounded(0);
        let sp = sh48_from_calibrate.0;
        thread::spawn(move || {
            let mut sh48 = OpticalPathToSH48::new();
            let __m = sh48.build().calibrate();
            sp.send(__m).unwrap();
        });
        // ---------------------------------------------------

        let rp = sh48_from_calibrate.1;
        //println!("{:?}", rp.recv());
        let __m = rp.recv().unwrap();

        //println!("{:?}", r.recv());

        let n_c: usize = (sh48.sensor.n_side_lenslet.pow(2) as usize) * 2 * sh48.sensor.n_sensor as usize;
        let gain = 0.5_f32;
        let mut _c_e: Array2<f32>;
        let mut rbm = Array2::<f32>::zeros((14, 6));
        let mut bm = Array2::<f32>::zeros((7, sh48.gmt.m1_n_mode as usize));

        loop {

            let gstate = match r.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            sh48.gmt.update(&gstate);

            sh48.propagate_src();
            sh48.sensor.process();
            let slopes = Array::from_shape_vec((n_c, 1), sh48.sensor.centroids.clone()).unwrap();

            _c_e = __m.dot(&slopes);
            let q = gain
                * _c_e
                .slice(s![..84, ..])
                .to_owned()
                .into_shape((14, 6))
                .unwrap();
            rbm -= &q.view();
            if sh48.gmt.m1_n_mode > 0 {
                let q = gain
                    * _c_e
                    .slice(s![84.., ..])
                    .to_owned()
                    .into_shape((7, sh48.gmt.m1_n_mode as usize))
                    .unwrap();
                bm -= &q.view();
            }

            let mut gstate = GmtState {
                rbm: Array2::<f32>::zeros((14, 6)),
                bm: Array2::<f32>::zeros((7, sh48.gmt.m1_n_mode as usize)),
            };
            gstate.rbm = rbm.clone();
            gstate.bm = bm.clone();
            s.send(gstate).unwrap();
 
        }

        let gstate = GmtState {
            rbm: Array2::<f32>::zeros((14, 6)),
            bm: Array2::<f32>::zeros((7, sh48.gmt.m1_n_mode as usize)),
        };
        drop(sh48);
        s.send(gstate).unwrap();
    });
    // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let r = plant_to_science.1;
    let s = science_plan_chat.0;
    thread::spawn(move || {
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

            let gstate = match r.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            gmt.update(&gstate);

            src_wfe_rms = src.through(&mut gmt).wfe_rms_10e(-9)[0];
            onaxis_pssn = src_pssn.reset(&mut src);

            println!(
                "@science => WFE RMS: {}nm; PSSn: {}",
                &src_wfe_rms.to_string(),
                &onaxis_pssn.to_string()
            );
            //thread::sleep(Duration::from_secs(1));
        }

        drop(gmt);
        drop(src);

        s.send("Bye from SCIENCE!".to_string()).unwrap();
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
    
    let s_science = plant_to_science.0;
    let s_sh48 = plant_to_sh48.0;
    let mut gstate = GmtState {
        rbm: Array2::<f32>::zeros((14, 6)),
        bm: Array2::<f32>::zeros((7, m1_n_mode as usize)),
    };
    gstate.rbm = gstate0.rbm.clone();
    gstate.bm = gstate0.bm.clone();

    let r = sh48_plan_chat.1;
     for _ in 0..20 {

        s_science.send(gstate.clone()).unwrap();
        s_sh48.send(gstate.clone()).unwrap();

        let new_gstate = match r.recv() {
            Ok(state) => state,
            Err(_) => break,
        };
        gstate.rbm = gstate0.rbm.clone() + new_gstate.rbm;
        gstate.bm = gstate0.bm.clone() + new_gstate.bm;

    }

    drop(s_science);
    drop(s_sh48);

    r.recv().unwrap();
    let r = science_plan_chat.1;
    println!("{:?}", r.recv());
}
