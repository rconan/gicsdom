use ceo::{CEO, element::ATMOSPHERE, Conversion};

fn main() {
    let n_px = 769;
    let atm_blueprint = CEO::<ATMOSPHERE>::new()
        .remove_turbulence_layer(0)
        .set_ray_tracing(
            25.5,
            n_px as i32,
            20f32.from_arcmin(),
            20f32,
            Some("glao_fiducial_atmosphere.bin".to_owned()),
            Some(20),
        )
        .build();
}
