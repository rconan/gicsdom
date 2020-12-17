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
        .point_of_origin(m1.height_at([10f64,0f64,0f64]))
    //        .polar_direction_vector(Arcsecond(1f64).to_radians(), 0f64)
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
}
