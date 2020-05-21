use csv::{Reader, Writer};
use gicsdom::astrotools;
use gicsdom::ceo::{set_gpu, Conversion, Gmt, PSSn, Source};
use rayon;
use rayon::prelude::*;
use serde::Deserialize;
use std::{f64, time::Instant};

#[derive(Debug, Deserialize, Default)]
struct ScienceField {
    zenith: f64,
    azimuth: f64,
}

#[derive(Debug, Deserialize, Default)]
struct StarField {
    index: usize,
    datetime: String,
    ALT: f64,
    AZ: f64,
    RA: f64,
    DEC: f64,
    Z1: f64,
    A1: f64,
    V1: f64,
    P1: f64,
    Z2: f64,
    A2: f64,
    V2: f64,
    P2: f64,
    Z3: f64,
    A3: f64,
    V3: f64,
    P3: f64,
    Z4: f64,
    A4: f64,
    V4: f64,
    P4: f64,
    r0: f64,
    L0: f64,
}

fn x_rotate(o: f64, w: Vec<f64>) -> Vec<f64> {
    let (s, c) = o.sin_cos();
    let (x, y, z) = (w[0], w[1], w[2]);
    vec![x, c * y - s * z, s * y + c * z]
}
fn y_rotate(o: f64, w: Vec<f64>) -> Vec<f64> {
    let (s, c) = o.sin_cos();
    let (x, y, z) = (w[0], w[1], w[2]);
    vec![c * x + s * z, y, -s * x + c * z]
}
fn z_rotate(o: f64, w: Vec<f64>) -> Vec<f64> {
    let (s, c) = o.sin_cos();
    let (x, y, z) = (w[0], w[1], w[2]);
    vec![c * x - s * y, s * x + c * y, z]
}
fn anti_gravity(alt: f64) -> Vec<f64> {
    let alpha = -13.522_f64.to_radians();
    let za = f64::consts::FRAC_PI_2 - alt;
    let mut ag: Vec<f64> = Vec::with_capacity(21);
    for k in 0..6 {
        let gamma = 60_f64.to_radians() * k as f64;
        let u = x_rotate(za, vec![0f64, 0f64, 1f64]);
        let v = z_rotate(gamma, u);
        let mut w = x_rotate(alpha, v);
        w[2] -= 1_f64;
        ag.append(&mut w);
    }
    let mut u = x_rotate(za, vec![0f64, 0f64, 1f64]);
    u[2] -= 1_f64;
    ag.append(&mut u);
    ag
}
#[test]
fn m1_local_gravity() {
    let mut ag = anti_gravity(90_f64.to_radians());
    let mut e = 0_f64;
    for k in 0..7 {
        let ag_k = &ag[3 * k..3 * (k + 1)];
        println!(
            "[{}]",
            ag_k.iter()
                .map(|x| format!("{:.9}", x))
                .collect::<Vec<String>>()
                .as_slice()
                .join(",")
        );
        if k < 6 {
            e += (ag_k[1] - 0.23381871).powi(2) + (ag_k[2] - -0.02771979).powi(2);
        }
    }
    println!("e: {}", e.sqrt());
    assert!(e < 1e-6);
}

fn pssn_temporal_stability(
    field: &StarField,
    sampling_time: f64,
) -> Result<(f32, f32, f32), String> {
    let mut science_field: Vec<ScienceField> = vec![];
    let mut rdr = Reader::from_path("KPP_field_sampler.csv").unwrap();
    for result in rdr.deserialize() {
        science_field.push(result.unwrap());
    }
    //println!("{:?}", science_field);
    let mut src = Source::new(science_field.len() as i32, 25.5, 401);
    let zen = science_field
        .iter()
        .map(|x| x.zenith.from_arcmin() as f32)
        .collect::<Vec<f32>>();
    let azi = science_field
        .iter()
        .map(|x| x.azimuth.to_radians() as f32)
        .collect::<Vec<f32>>();
    let mag = vec![0f32; science_field.len()];
    src.build("Vs", zen, azi, mag);

    let mut gmt = Gmt::new();
    gmt.build_m1("M1_printthrough", 3).build_m2(None);
    src.through(&mut gmt).xpupil();
    //println!("WFE RMS: {:.3}nm", src.wfe_rms_10e(-9)[0]);

    let mut pssn = PSSn::new(0.15, 25.0);
    pssn.build(&mut src);

    //println!("{}", field.datetime);
    let star_field_time = astrotools::Time::from_iso_8601(&field.datetime);
    //println!("{:?}", star_field_time);

    let target = astrotools::SkyCoordinates::new(field.RA, field.DEC);

    //let sampling_time = 60_f64;
    let duration = 8_f64 * 15_f64 * 60_f64;
    let peek_time = 15_f64 * 60_f64;
    let peek_lag = (peek_time / sampling_time).round() as usize;
    //println!("Peek lag: {}", peek_lag);
    let mut obs = astrotools::Observation::new(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        star_field_time,
        target,
        sampling_time,
        duration,
    );
    let (alt, az) = obs.altaz();
    /*
    println!(
        " == {} - [{:.3},{:.3}]",
        obs.utc.datetime(),
        alt.to_degrees(),
        az.to_degrees()
    );
    */
    pssn.reset(&mut src);
    let psu = pssn.spatial_uniformity();
    //println!("PSSn: {} - PSU: {:.3}%", pssn, psu);

    let mut sorted_field_pssn: Vec<Vec<f32>> = vec![Vec::with_capacity(8); 21];
    let mut psu: Vec<f32> = Vec::with_capacity(8);

    while let Some(_) = obs.next() {
        let (alt, az) = obs.altaz();
        let alt_deg = alt.to_degrees();
        if alt_deg > 89.5 || alt_deg < 30.0 {
            return Err("Elevation out of allowed range [30,89.5]!".to_string());
        }
        let mut ag = anti_gravity(alt);
        gmt.set_m1_modes(&mut ag);
        src.through(&mut gmt).xpupil();
        pssn.accumulate(&mut src);
        if (obs.step % peek_lag) == 0 {
            /*
            println!(
                " == {} - [{:.3},{:.3}]",
                obs.utc.datetime(),
                alt_deg,
                az.to_degrees()
            );
            */
            pssn.peek(&mut src);
            psu.push(pssn.spatial_uniformity());
            //println!("PSSn: {} - PSU: {:.3}%", pssn, psu.last().unwrap());

            for k in 0..21 {
                sorted_field_pssn[k].push(pssn.estimates[k])
            }
        }
    }

    drop(src);
    drop(gmt);
    drop(pssn);

    let pssn_on_axis = sorted_field_pssn[0][0];
    let mut median_field_pssn = vec![];
    for k in 0..21 {
        sorted_field_pssn[k].sort_by(|a, b| a.partial_cmp(b).unwrap());
        median_field_pssn.push(0.5 * (sorted_field_pssn[k][3] + sorted_field_pssn[k][4]));
    }
    median_field_pssn.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean_median_field_pssn = median_field_pssn.iter().sum::<f32>() / 21 as f32;
    let pts = 100. * (median_field_pssn.last().unwrap() - median_field_pssn.first().unwrap())
        / mean_median_field_pssn;
    //println!("{:?}", median_field_pssn);
    //println!("{:?}", mean_median_field_pssn);
    //println!("PTS: {:.3}%", pts);

    Ok((pssn_on_axis, psu[0], pts))
}
fn main() {
    let mut star_field: Vec<StarField> = vec![];
    let mut rdr = Reader::from_path("guide_stars05.csv").expect("File with star fields missing!");
    for result in rdr.deserialize() {
        match result {
            Ok(record) => star_field.push(record),
            Err(err) => println!("{}", err),
        }
    }
    println!("# star field: {}", star_field.len());

    let now = Instant::now();

    let sampling_time = 10_f64;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(32)
        .build()
        .unwrap();
    //let star_field_id = 0;
  
    //star_field.truncate(8);
    let results = pool.install(|| {
        star_field
            .par_iter()
            .map(|field| {
                let thread_id = pool.current_thread_index().unwrap();
                //        println!("Thread ID#{}",thread_id);
                set_gpu((thread_id % 4) as i32);
                pssn_temporal_stability(&field, sampling_time)
            })
            .collect::<Vec<Result<(f32, f32, f32), String>>>()
    });
//    println!("{:?}", results);

    println!("ET: {}s", now.elapsed().as_secs());

    let mut wtr = Writer::from_path("KPP_M1_Printthrough.csv").unwrap();
    wtr.write_record(&["index","oa","psu","pts"]); 
    for (r,f) in results.iter().zip(star_field) {
        match r {
            Ok((pssn, psu, pts)) => wtr.serialize((f.index,pssn,psu,pts)).unwrap(),
            Err(e) => println!("Star field #{} out of elevation range!",f.index),
        }
    }
    wtr.flush().unwrap();
}
