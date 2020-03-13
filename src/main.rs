//#[macro_use]
extern crate crossbeam_channel;

use console::Term;
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use csv::Reader;
use gicsdom::ceo;
use gicsdom::ceo::GmtState;
use gicsdom::{Observation, OpticalPathToDSH48, OpticalPathToSH48, SkyCoordinates};
//use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use ndarray::{s, stack, Array, Array2, Axis};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use ndarray_npy::{NpzReader, NpzWriter};
use serde::Deserialize;
use std::fs::File;
use std::io::ErrorKind;
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

#[derive(Debug, Deserialize)]
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
#[derive(Debug, Deserialize)]
struct ScienceField {
    zenith: f64,
    azimuth: f64,
}

fn main() {
    let chat_plant_sh48: [Session<GmtState, Vec<f32>>; 3] = [
        Session::new(bounded(0), bounded(0)),
        Session::new(bounded(0), bounded(0)),
        Session::new(bounded(0), bounded(0)),
    ];
    let chat_plant_science: Session<GmtState, &str> = Session::new(bounded(0), unbounded());

    let mut field = StarField {
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
    let probe_ids = [field.p2 as i32, field.p3 as i32, field.p4 as i32];
    let probe_sh48_gs_mag = [field.v2, field.v3, field.v4];

    let sampling = 30.0f64;

    let term = Term::buffered_stdout();

    // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let mut h_48 = vec![];
    for (i, chat_plant_sh48_) in chat_plant_sh48.iter().enumerate() {
        let sh48_datetime = datetime.clone();
        let sh48_term = term.clone();
        let (to_plant, from_plant) = chat_plant_sh48_.dual_channel();
        h_48.push(thread::spawn(move || {
            ceo::set_gpu(i as i32+1);

            let mut atm = ceo::Atmosphere::new();
            atm.load_from_json("/home/ubuntu/DATA/gmtAtmosphereL025_1579821046.json").unwrap();
            let atm_sampling = 1e-1f64;

            let exposure_time = 30.0;
            let readout_noise_rms = 0.5;
            let n_background_photon = 0.0;
            let noise_factor = 1.4142;

            let mut sh48 =
                OpticalPathToDSH48::new(&sh48_datetime, telescope, probes[i], probe_ids[i]);
            sh48.build(8, Some(24), None, probe_sh48_gs_mag[i]);

            loop {
                let gstate = match from_plant.recv() {
                    Ok(state) => state,
                    Err(_) => break,
                };
                sh48.gmt.update(&gstate);

                let mut stopwatch = 0.0f64;
                while stopwatch<sampling {
                    let (z, a, da, p) = sh48.update(atm_sampling, &mut atm).local();
                    stopwatch += atm_sampling;
                    atm.secs += atm_sampling;
  
                    sh48_term.move_cursor_to(0,46+i).unwrap();
                    sh48_term
                        .write_line(&format!(
                            "SH48 #{}  {:>4.2}' {:>+7.2} {:>+6.2} {:>+6.2} [{:5.2}]",
                            sh48.probe_id, z, a,da,p,
                            stopwatch
                        ))
                        .unwrap();
                    sh48_term.flush().unwrap();
  
                }

                sh48.sensor
                    .readout(
                        exposure_time,
                        readout_noise_rms,
                        n_background_photon,
                        noise_factor,
                    )
                    .process();

                to_plant.send(sh48.sensor.centroids.clone()).unwrap();
            }
            drop(sh48);
        }));
    }
        // SH48 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // SCIENCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    let mut science_field = ScienceField {
        zenith: 0.0,
        azimuth: 0.0,
    };
    let mut zen = vec![];
    let mut azi = vec![];
    let mut rdr = Reader::from_path("KPP_field_sampler.csv").unwrap();
    for result in rdr.deserialize() {
        science_field = result.unwrap();
        zen.push((science_field.zenith.to_radians()/60.) as f32);
        azi.push(science_field.azimuth.to_radians() as f32);
    }

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
        let n_src = zen.len();
        let mut src = ceo::Source::new(n_src as i32, pupil_size, pupil_sampling);
        src.build("V", zen, azi, vec![0.0;n_src]);

        let mut src_wfe_rms = src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0];
        let mut onaxis_pssn: f32 = 0.0;
        let mut pssn_spatial_uniformity: f32 = 0.0;
        let mut src_pssn = ceo::PSSn::new(15e-2, 25.0, 0.0);
        src_pssn.build(&mut src);

        loop {
            let gstate = match from_plant.recv() {
                Ok(state) => state,
                Err(_) => break,
            };

            gmt.update(&gstate);

            src_wfe_rms = src.through(&mut gmt).xpupil().wfe_rms_10e(-9)[0];
            src_pssn.peek(&mut src);
            onaxis_pssn = src_pssn.estimates[0];
            pssn_spatial_uniformity = src_pssn.spatial_uniformity();

            sci_term.move_cursor_to(0,45).unwrap();
            sci_term
                .write_line(&format!(
                    "WFE RMS: {:.1} nm; PSSn: {:>6.4} {:>6.4}%",
                    src_wfe_rms, onaxis_pssn, pssn_spatial_uniformity
                ))
                .unwrap();
            sci_term.flush().unwrap();
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
/*    gstate0.rbm[(0, 0)] = 1e-5;
    gstate0.rbm[(8, 4)] = 1e-6;
    gstate0.rbm[(2, 3)] = -1e-6;
    gstate0.rbm[(11, 1)] = -1e-5;
    gstate0.bm[[0, 0]] = 1e-5;
    gstate0.bm[[1, 0]] = 1e-5;
    gstate0.bm[[2, 1]] = 1e-5;
*/
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


    let f = File::open("recon.npz").unwrap_or_else(|error| {
        if error.kind() == ErrorKind::NotFound {
            File::create("recon.npz").unwrap_or_else(|error| {
                panic!("Problem creating the file: {:?}", error);
            })
        } else {
            panic!("Problem opening the file: {:?}", error);
        }
    });

    let mut npz = NpzReader::new(f).unwrap();
    let __m: Array2<f32> = npz.by_name("m").unwrap();

/*
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

    let mut npz = NpzWriter::new(File::create("recon.npz").unwrap());
    npz.add_array("m", &__m).unwrap();
*/
    //    sp.send(__m).unwrap();
    //});
    //let rp = sh48_from_calibrate.1;
    //println!("{:?}", rp.recv());
    //let __m = rp.recv().unwrap();
    // ---------------------------------------------------
    thread::sleep(Duration::from_secs(20));

    let mut obs = Observation::new(&datetime, -29.04, -70.682);
    let pointing = SkyCoordinates::new(telescope.0, telescope.1);
    let (alt0, az0, pa0) = pointing.alt_az_parallactic_deg(&obs);

    term.move_cursor_to(0,42).unwrap();
    term.write_line(&format!("{:=>75}", "",)).unwrap();
    let n = (15. * 60. / sampling).ceil() as u64;
    let mut cum_et = 0u64;
    let now = Instant::now();
    for k in 0..n {
        let now = Instant::now();

        let (alt, az, pa) = pointing.alt_az_parallactic_deg(&obs);
        term.move_cursor_to(0,43).unwrap();
        term.write_line(&format!(
            "  {} ({:+>6.2},{:+>5.2}) ({:+>6.2},{:+>5.2}) ({:+>6.2},{:+>5.2}) [{}/{}]",
            obs.datetime, alt, alt-alt0, az, az-az0, pa, pa-pa0, k+1,n
        ))
        .unwrap();
        term.write_line(&format!("{:->75}", "")).unwrap();
        term.flush().unwrap();

        to_science.send(gstate.clone()).unwrap();

        let mut q: Vec<f32> = Vec::with_capacity((n_c) as usize);


        for chat_plant_sh48_ in chat_plant_sh48.iter() {
            let (to_sh48, from_sh48) = chat_plant_sh48_.main_channel();
            to_sh48.send(gstate.clone()).unwrap();
        }
        for chat_plant_sh48_ in chat_plant_sh48.iter() {
            let (to_sh48, from_sh48) = chat_plant_sh48_.main_channel();

            let mut centroids = match from_sh48.recv() {
                Ok(state) => state,
                Err(_) => break,
            };
            q.append(&mut centroids);
        }
        let slopes = Array::from_shape_vec((n_c, 1), q).unwrap();
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

        term.move_cursor_to(0,49).unwrap();
        term.write_line(&format!("{:->75}", "")).unwrap();
        let et = now.elapsed().as_secs();
        cum_et += et;
        let eta = (n-k-1)*cum_et/(k+1);
        term.write_line(&format!(" Step: {:>4}s ; Total: {:>5}s ; ETA: {:>5}s", et, cum_et, eta))
            .unwrap();
        term.write_line(&format!("{:=>75}", "",)).unwrap();
    }
//    let et = now.elapsed().as_secs();
//    println!("ET: {}s",et);
/*
    term.move_cursor_down(6).unwrap();
    term.write_line(&format!("{:->60}", "")).unwrap();
    let et = now.elapsed().as_secs();
    let sim_rate = n / et;
    term.write_line(&format!(" Done in {}s ({:.2}steps/s)", et, sim_rate))
        .unwrap();
    term.write_line(&format!("{:=>60}", "",)).unwrap();
    term.flush().unwrap();
*/
    drop(to_science);
    for chat_plant_sh48_ in &chat_plant_sh48 {
        let (to_plant, from_plant) = chat_plant_sh48_.dual_channel();
        drop(from_plant);
        drop(to_plant);
        let (to_sh48, from_sh48) = chat_plant_sh48_.main_channel();
        drop(to_sh48);
            drop(from_sh48);
    }
    for h in h_48 {
        h.join().unwrap();
    }

    h_science.join().unwrap();
}
