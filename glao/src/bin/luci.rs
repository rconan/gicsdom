use ceo::Conversion;
//use glao::glao_sys;
use glao::glao_sys::{GlaoSys, ScienceField};
use glao::system::Cn2;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle};
use log::LevelFilter;
use rayon::prelude::*;
use simple_logger::SimpleLogger;
use std::f32;

fn main() {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .init()
        .unwrap();

    let mut cn2_reader = csv::Reader::from_path("glao_cn2.csv").unwrap();
    let mut cn2_profiles: Vec<Cn2> = vec![];
    for result in cn2_reader.deserialize() {
        cn2_profiles.push(result.unwrap());
    }
    let secz = 1f32 / 30f32.to_radians().cos();

    let m1_polishing_error = false;
    let n_sample: usize = 1000;

    let n_gpu = 8 as usize;
    let n_thread = n_gpu;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_thread)
        .build()
        .unwrap();


    pool.install(|| {
        let mpb = MultiProgress::new();
        let pb_main = mpb.add(ProgressBar::new(cn2_profiles.len() as u64));
        pb_main.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.green/yellow} {pos:>4}/{len:4}"),
        );
        cn2_profiles
            .par_iter()
            .progress_with(pb_main)
            .for_each(|cn2_prof| {
                let thread_id = pool.current_thread_index().unwrap();
                ceo::set_gpu((thread_id % n_gpu) as i32);

                log::info!("LUCI CN2 #{}", cn2_prof.idx);

                let turb_cn2_height = [40f32, 125f32, 350f32, 1500f32, 4000f32, 8000f32, 16000f32]
                    .iter()
                    .map(|x| x * secz)
                    .collect::<Vec<f32>>();
                let turb_cn2_xi0 = [
                    cn2_prof.m40,
                    cn2_prof.m125,
                    cn2_prof.m350,
                    cn2_prof.m1500,
                    cn2_prof.m4000,
                    cn2_prof.m8000,
                    cn2_prof.m16000,
                ];

                let atm_r0 = 0.9759 * 500e-9 / cn2_prof.dimm.from_arcsec();
                let atm_oscale = 25f32;
                log::info!("Atmosphere: r0={} , L0={}", atm_r0, atm_oscale);

                let n_px = (4.0 * 25.5
                    / ceo::PSSn::<ceo::pssn::TelescopeError>::r0_at_z(
                        atm_r0 as f32,
                        30f32.to_radians(),
                    ))
                .ceil()
                .min(1024f32) as usize;
                log::info!("Source pupil sampling: {}", n_px);

                let mut atm = ceo::Atmosphere::new();
                let mut science =
                    ScienceField::delaunay_21("Vs", n_px, Some((atm_r0 as f32, atm_oscale as f32)));
                science.build();
                atm.build(
                    science.pssn.r0(),
                    science.pssn.oscale,
                    7,
                    turb_cn2_height,
                    turb_cn2_xi0.to_vec(),
                    vec![0f32; 7],
                    vec![0f32; 7],
                );
                let mut glao_4gs = GlaoSys::default(&mut atm, &mut science);
                glao_4gs.build(6f32.from_arcmin(), 70, 0.5).calibration();
                if m1_polishing_error {
                    glao_4gs.s12 = Some(((2, 7), (1, 1)));
                    //glao_sys::m1_polishing_wavefront_error(&mut glao_4gs);
                }
                let pb = mpb.add(ProgressBar::new(n_sample as u64));
                pb.set_style(ProgressStyle::default_bar().template(&format!(
                    "{:>6}: {}",
                    cn2_prof.idx, "[{eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}"
                )));
                let mut pssn_previous = 0f32;
                let mut tol_count = 0;
                let mut tol = 0f32;
                for k in 0..n_sample {
                    pb.inc(1);
                    glao_4gs.next();
                    if k % 25 == 0 {
                        let pssn_current = glao_4gs.peek()[12];
                        tol = (pssn_previous - pssn_current).abs();
                        pb.set_message(&format!("Tol.: {:>7.4}", tol));
                        if tol < 1e-4 {
                            tol_count += 1;
                        } else {
                            tol_count = 0;
                        }
                        pssn_previous = pssn_current;
                    }
                    if tol_count == 2 {
                        science.pssn_nsample_tol = Some((k, tol));
                        break;
                    }
                }
                pb.finish();
                science.wrap_up();
                if science.pssn_nsample_tol.is_none() {
                    let pssn_current = science.pssn.estimates[12];
                    tol = (pssn_previous - pssn_current).abs();
                    science.pssn_nsample_tol = Some((n_sample, tol));
                }
                let filename = format!(
                    "Results/glao_luci_{:04}cn2_{:04}.pkl",
                    cn2_prof.idx, n_sample
                );
                science.dump(&filename).unwrap();
            });
        mpb.join_and_clear().unwrap();
    });

}
