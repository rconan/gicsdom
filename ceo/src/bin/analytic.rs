use ceo::analytic::*;
use ceo::SkyAngle::*;

fn main() {
    let m1 = Conic::gmt_m1();
    let m2 = Conic::gmt_m2();

    println!("CHIEF RAY:");
    let mut ray = new_ray()
        .polar_direction_vector(Arcminute(1f64).to_radians(), 0f64)
        .build();
    println!("Init   : {}", ray);
    m1.reflect(&mut ray);
    println!("Reflect: {}", ray);
    ray.trace_to(&m2);
    println!("Trace  : {}", ray);
    m2.reflect(&mut ray);
    println!("Reflect: {}", ray);
    let z_xpupil = ray.solve_for_z(0f64, 0f64);
    println!("Exit pupil: {:.9}m", z_xpupil);

    println!("MARGINAL RAY:");
    let mut ray = new_ray()
        .point_of_origin(m1.height_at([10f64, 0f64, 0f64]))
        //        .polar_direction_vector(Arcminute(10f64).to_radians(), 0f64)
        .build();
    println!("Init   : {}", ray);
    m1.reflect(&mut ray);
    println!("Reflect: {}", ray);
    let z_gregorian = ray.solve_for_z(0f64, 0f64);
    println!("Gregorian focus: {:.9}m", z_gregorian);
    ray.trace_to(&m2);
    println!("Trace  : {}", ray);
    m2.reflect(&mut ray);
    println!("Reflect: {}", ray);
    let z_focal_plane = ray.solve_for_z(0f64, 0f64);
    println!("Focal plane: {:.9}m", z_focal_plane);

    let gmt = Gmt::new();
    let m: Vec<Vector> = (1..11).map(|x| [x as f64, 0., 0.]).collect();
    let rays = gmt.focal_point(m.clone(), Arcminute(10f64).to_radians(), 0f64);
    rays.iter()
        .enumerate()
        .for_each(|x| println!("#{:3}: {}", x.0, x.1));

    for z in 1..11 {
        let rays = gmt.focal_point(m.clone(), Arcminute(z as f64).to_radians(), 0f64);
        let z_mean = rays.iter().map(|r| r.p[2]).sum::<f64>() / rays.len() as f64;
        let z_rms = (rays.iter().map(|r| (r.p[2] - z_mean).powf(2.)).sum::<f64>()
            / rays.len() as f64)
            .sqrt();
        let mut p = rays[0].p;
        p[2] = z_mean - z_focal_plane;
        let r = -0.5 * p.norm_square() / p[2];
        println!(
            "z: {:2}': {:.9}m +/- {:3.0}nm ; curvature: {:.9}m",
            z,
            z_mean,
            z_rms * 1e9,
            r
        );
    }
}
