#![allow(non_snake_case)]
use csv::{Reader, Writer};
use gicsdom::astrotools;
use rayon;
use rayon::prelude::*;
use serde::Deserialize;
use std::{f64, time::Instant};


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

    let field = &star_field[100];

    let star_field_time = astrotools::Time::from_iso_8601(&field.datetime);
    //println!("{:?}", star_field_time);

    let target = astrotools::SkyCoordinates::new(field.RA, field.DEC);

    let sampling_time = 1_f64;
    let duration = 10_f64;
    //println!("Peek lag: {}", peek_lag);
    let mut obs = astrotools::Observation::new(
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
        star_field_time,
        target,
        sampling_time,
        duration,
    );

    let mut wtr = Writer::from_path("mount_trajectory.csv").unwrap();
    wtr.write_record(&["alt", "az"]).unwrap();

    while let Some(_) = obs.next() {
        let (alt, az) = obs.altaz_deg();
        println!("alt: {} ; az: {}",alt,az);
        wtr.serialize((alt,az)).unwrap()
    }
    wtr.flush().unwrap();
}
