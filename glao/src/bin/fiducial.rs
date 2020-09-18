use ceo::Conversion;
use glao::glao_sys;
use glao::glao_sys::{GlaoSys, ScienceField};
use indicatif::{ProgressBar, ProgressStyle};
use log::LevelFilter;
use simple_logger::SimpleLogger;

fn main() {
    SimpleLogger::new()
        .with_level(LevelFilter::Info)
        .init()
        .unwrap();

    let m1_polishing_error = false;

    let mut atm = ceo::Atmosphere::new();
    let mut science = ScienceField::on_axis("Vs", 1024, None);
    science.build();
    atm.gmt_build(science.pssn.r0(), science.pssn.oscale);
    let mut glao_4gs = GlaoSys::default(&mut atm, &mut science);
    glao_4gs.build(6f32.from_arcmin(), 70, 0.5).calibration();
    if m1_polishing_error {
        glao_4gs.s12 = Some(((2, 7), (1, 1)));
        glao_sys::m1_polishing_wavefront_error(&mut glao_4gs);
    }

    let n_sample = 1;
    let pb = ProgressBar::new(n_sample as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7}"),
    );
    let (atm_pssn, pssn, atm_fwhm_x, atm_fwhm, glao_fwhm) = {
        for k in &mut glao_4gs {
            pb.inc(1);
            if k == n_sample {
                break;
            }
        }
        pb.finish();
        science.wrap_up().results()
    };
    println!(
        "PSSN: {:?}/{:?}={} ; FWHM [arcsec]: {:.3}/{:.3}/{:.3}",
        pssn,
        atm_pssn,
        pssn[0] / atm_pssn[0],
        atm_fwhm_x,
        atm_fwhm[0],
        glao_fwhm[0]
    );
}
