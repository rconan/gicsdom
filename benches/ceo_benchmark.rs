extern crate gicsdom;
use gicsdom::ceo;
use crate::gicsdom::ceo::Propagation;
use std::f32;
use std::time::Duration;

#[macro_use]
extern crate criterion;

use criterion::Criterion;

fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 1,
        1 => 1,
        n => fibonacci(n-1) + fibonacci(n-2),
    }
}


fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib 20", |b| b.iter(|| fibonacci(20)));
}

fn ceo_atmosphere_get_screen(c: &mut Criterion) {
    let mut gmt = ceo::Gmt::new(0,None);
    gmt.build();
    let mut atm = ceo::Atmosphere::new();
    atm.build(0.15,30.0,26.0,401,0.0,1.0,"/home/ubuntu/DATA/test.bin",1);

    let mut group = c.benchmark_group("CEO ATM");
    group.sample_size(10);
    for k in 1..11{
        let n = k*100+1;
        let mut src = ceo::Source::new(1,25.5,n);
        src.build("V",vec![0.0f32],vec![0.0f32],vec![0.0f32]);
        let fname = format!("n={}",n);
        group.bench_function(fname, |b| b.iter(|| {
            src.through(&mut gmt);
            atm.propagate(&mut src);
        })
        );
    }
}

fn ceo_shack_hartmann(c: &mut Criterion) {
    let mut gmt = ceo::Gmt::new(0,None);
    gmt.build();
    let mut wfs = ceo::ShackHartmann::new(1,48,16,25.5/48.0);
    wfs.build(8,Some(24),None);
    let mut src = wfs.new_guide_stars();
    src.build("V",vec![0.0f32],vec![0.0f32],vec![0.0f32]);
//    let mut atm = ceo::Atmosphere::new();
//    atm.build(0.15,30.0,26.0,401,0.0,1.0,"/home/ubuntu/DATA/test.bin",1);

    let mut group = c.benchmark_group("CEO SH WFS");
    group.sample_size(10);
    group.bench_function("SH48", |b| b.iter(|| {
        src.through(&mut gmt).through(&mut wfs);
    })
    );
}

fn ceo_shack_hartmann_with_atmosphere(c: &mut Criterion) {
    let mut gmt = ceo::Gmt::new(0,None);
    gmt.build();
    let mut wfs = ceo::ShackHartmann::new(1,48,16,25.5/48.0);
    wfs.build(8,Some(24),None);
    let mut src = wfs.new_guide_stars();
    src.build("V",vec![0.0f32],vec![0.0f32],vec![0.0f32]);
    let mut atm = ceo::Atmosphere::new();
    atm.build(            0.15,
            30.0,
            26.0,
            1300,
            20.0 * f32::consts::PI / 180. / 60.,
            15.0,
            "/home/ubuntu/DATA/gmtAtmosphereL025_1579821046.bin",
            80,
    );

    let mut group = c.benchmark_group("CEO SH WFS ATM");
    group.sample_size(10);
    group.bench_function("SH48_ATM", |b| b.iter(|| {
        src.through(&mut gmt);
        src.is_opd_to_phase = false;
        atm.propagate(&mut src);
        src.is_opd_to_phase = true;
        src.through(&mut wfs);
    })
    );
}

fn ceo_shack_hartmann_with_atmosphere_integrated(c: &mut Criterion) {
    let mut gmt = ceo::Gmt::new(0,None);
    gmt.build();
    let mut wfs = ceo::ShackHartmann::new(1,48,16,25.5/48.0);
    wfs.build(8,Some(24),None);
    let mut src = wfs.new_guide_stars();
    src.build("V",vec![0.0f32],vec![0.0f32],vec![0.0f32]);
    let mut atm = ceo::Atmosphere::new();
    atm.build(            0.15,
            30.0,
            26.0,
            1300,
            20.0 * f32::consts::PI / 180. / 60.,
            15.0,
            "/home/ubuntu/DATA/gmtAtmosphereL025_1579821046.bin",
            80,
    );

    let mut group = c.benchmark_group("CEO SH WFS ATM INT.");
    group.sample_size(10);

    group.bench_function("SH48_ATM", |b| b.iter(|| {
        src.through(&mut gmt);
        src.is_opd_to_phase = false;
        atm.propagate(&mut src);
        src.is_opd_to_phase = true;
        src.through(&mut wfs);
    })
    );

    group.measurement_time(Duration::from_secs(600));

    group.bench_function("SH48_INT_ATM", |b| b.iter(|| {
        src.through(&mut gmt);
        atm.secs = 0.0;
        let atm_sampling = 0.01;
        let inc_secs = 30.0;
        let n = (inc_secs/atm_sampling) as u32;
        for k in 0..n {
            src.is_opd_to_phase = false;
            atm.propagate(&mut src);
            src.is_opd_to_phase = true;
            src.through(&mut wfs);
            atm.secs += atm_sampling;
        }
    })
    );
}


criterion_group!(benches,
                 ceo_shack_hartmann_with_atmosphere_integrated);
criterion_main!(benches);
