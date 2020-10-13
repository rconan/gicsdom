use ceo::{Atmosphere,Conversion};
use glao::glao_sys::ScienceField;

fn main() {
    let n_px = 769;
    let secz = 1f32 / 30f32.to_radians().cos();
    let turb_cn2_height = [275, 425, 1250, 4000, 8000, 13000]
        .iter()
        .map(|x| *x as f32 * secz)
        .collect::<Vec<f32>>();
    let turb_cn2_xi0 = [0.0874, 0.0666, 0.3498, 0.2273, 0.0681, 0.0751];
    let wind_speed = [5.7964, 5.8942, 6.6370, 13.2925, 34.8250, 29.4187];
    let wind_direction = [0.1441, 0.2177, 0.5672, 1.2584, 1.6266, 1.7462];
    let mut atm = Atmosphere::new();
    //let mut science =
    //    ScienceField::delaunay_21("Vs", n_px, None);
    let mut science = ScienceField::on_axis("Vs", n_px, None);
    science.build();
    atm.raytrace_build(
        science.pssn.r0(),
        science.pssn.oscale,
        6,
        turb_cn2_height,
        turb_cn2_xi0.to_vec(),
        wind_speed.to_vec(),
        wind_direction.to_vec(),
        25.5,
        n_px as i32,
        20f32.from_arcmin(),
        20f32,
        Some("glao_fiducial_atmosphere.bin"),
        Some(20)
    );
}
