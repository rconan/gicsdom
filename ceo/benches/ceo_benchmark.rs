use ceo::{
    element::*, shackhartmann::Geometric, CEOInit, Conversion, Gmt, ShackHartmann, Source, CEO,
};
use criterion::BenchmarkId;
use criterion::Throughput;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

#[inline]
fn sh48_ray_tracing(gmt: &mut Gmt, gs: &mut Source) {
    gs.through(gmt).xpupil();
}

#[inline]
fn sh48_wavefront_sensing(gmt: &mut Gmt, gs: &mut Source, wfs: &mut ShackHartmann<Geometric>) {
    gs.through(gmt).xpupil().through(wfs);
}

pub fn sh48_benchmark(c: &mut Criterion) {
    let mut gmt = CEO::<GMT>::new().build();
    let mut gs = (1..=4)
        .map(|n_sensor| {
            CEO::<SH48>::new()
                .set_n_sensor(n_sensor)
                .guide_stars()
                .set_on_ring(6f32.from_arcmin())
                .build()
        })
        .collect::<Vec<Source>>();
    let mut wfs = (1..=4)
        .map(|n_sensor| {
            CEO::<SH48>::new()
                .set_n_sensor(n_sensor)
                .build::<Geometric>()
        })
        .collect::<Vec<ShackHartmann<Geometric>>>();
    let mut group = c.benchmark_group("sh48_benchmark");
    for n_sensor in 1..=4 {
        //group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(
            BenchmarkId::new("Ray tracing", n_sensor),
            &n_sensor,
            |b, &n_sensor| b.iter(|| sh48_ray_tracing(&mut gmt, &mut gs[n_sensor - 1])),
        );
        group.bench_with_input(
            BenchmarkId::new("Wavefront sensing", n_sensor),
            &n_sensor,
            |b, &n_sensor| {
                b.iter(|| {
                    sh48_wavefront_sensing(&mut gmt, &mut gs[n_sensor - 1], &mut wfs[n_sensor - 1])
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, sh48_benchmark);
criterion_main!(benches);
