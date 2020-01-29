use gicsdom::optics_path_gsh48::OpticsPathGSH48;
use gicsdom::optics_path_sh48::OpticsPathSH48;

fn main() {
    let mut sh48 = OpticsPathSH48::new();
    sh48.cfg();
    sh48.goxp();
    println!("WFE RMS: {}nm", sh48.wfe_rms * 1e9);

    let mut gsh48 = OpticsPathGSH48::new();
    gsh48.cfg();
    let t_xyz: [f64; 3] = [1e-6, 0.0, 0.0];
    let r_xyz: [f64; 3] = [0.0, 1e-6, 0.0];
    for seg_id in 0..7 {
        gsh48.m1_update(t_xyz, r_xyz, seg_id);
    }
    gsh48.m1_status();
}
