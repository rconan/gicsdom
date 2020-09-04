use glao::system::System;

fn main() {
    ceo::set_gpu(1);
    let pupil_size = 25.5;
    let n_lenslet = 48;
    //let n_actuator = n_lenslet + 1;
    let n_px_lenslet = 16;
    let wfs_intensity_threshold = 0.5;

    let mut a_wfe_rms: Vec<f32> = vec![];
    for n_kl in (20..201).step_by(10) {
        println!("N_KL: {}", n_kl);

        let mut on_axis_sys = System::new(pupil_size, 1, n_lenslet, n_px_lenslet);
        on_axis_sys
            .gmt_build("bending modes", 27, n_kl)
            .wfs_build("V", vec![0f32], vec![0f32], vec![0f32])
            .wfs_calibrate(wfs_intensity_threshold)
            .through();

        let mut v = vec![0usize; (n_lenslet * n_lenslet) as usize];
        on_axis_sys
            .wfs
            .lenset_mask()
            .from_dev()
            .chunks((n_lenslet * n_lenslet) as usize)
            .for_each(|x| {
                /*let s = x
                    .iter()
                    .map(|x| if *x > 0f32 { 1 } else { 0 })
                    .sum::<usize>();
                println!("lenslet mask: {}", s);*/
                for k in 0..v.len() {
                    v[k] += if x[k] > 0f32 { 1usize } else { 0usize };
                }
            });
        let mask = v
            .iter()
            .map(|&x| if x > 0usize { 1u8 } else { 0u8 })
            .collect::<Vec<u8>>();
        let nnz = mask.iter().cloned().map(|x| x as usize).sum::<usize>();

        let mut full_calib = on_axis_sys.m2_mode_calibrate();

        let n = 7 * (n_kl - 1);
        let m = (n_lenslet * n_lenslet) as usize;
        let mut reduced_calib = full_calib
            .from_dev()
            .chunks(m as usize)
            .map(|x| {
                let mut r = vec![];
                for k in 0..m {
                    if mask[k] > 0 {
                        r.push(x[k])
                    }
                }
                r
            })
            .flatten()
            .collect::<Vec<f32>>();
        //println!("Reduced calibration: {}",reduced_calib.len());

        let mut calib: ceo::Cu<f32> = ceo::Cu::array(2 * nnz, n);
        calib.to_dev(&mut reduced_calib);
        calib.qr();

        let mut atm = ceo::Atmosphere::new();
        {
            let mut pssn: ceo::PSSn<ceo::pssn::AtmosphereTelescopeError> = ceo::PSSn::new();
            pssn.build(&mut on_axis_sys.gs);
            atm.gmt_build(pssn.r0(), pssn.oscale);
        }

        let n_sample = 10;
        let mut a_wfe_var = 0f32;
        for _ in 0..n_sample {
            on_axis_sys.gmt.reset();
            on_axis_sys.wfs.reset();
            on_axis_sys
                .gs
                .through(&mut on_axis_sys.gmt)
                .xpupil()
                .through(&mut atm)
                .through(&mut on_axis_sys.wfs);
            let wfe_rms_0 = on_axis_sys.gs.wfe_rms_10e(-9)[0];

            let mut red_s = on_axis_sys
                .wfs
                .centroids
                .from_dev()
                .chunks(m as usize)
                .map(|x| {
                    let mut r = vec![];
                    for k in 0..m {
                        if mask[k] > 0 {
                            r.push(x[k])
                        }
                    }
                    r
                })
                .flatten()
                .collect::<Vec<f32>>();

            let mut c: ceo::Cu<f32> = ceo::Cu::vector(red_s.len());
            c.to_dev(&mut red_s);
            let mut x = calib.qr_solve(&mut c);
            let h_x = x.from_dev();
            let mut kl_coefs = vec![vec![0f64; n_kl]; 7];
            let mut k = 0;
            for s in 0..7 {
                for a in 1..n_kl {
                    kl_coefs[s][a] -= h_x[k] as f64;
                    k += 1;
                }
            }
            let mut m = kl_coefs.clone().into_iter().flatten().collect::<Vec<f64>>();
            on_axis_sys.gmt.set_m2_modes(&mut m);
            on_axis_sys.gs.reset();
            on_axis_sys
                .gs
                .through(&mut on_axis_sys.gmt)
                .xpupil()
                .through(&mut atm);
            let wfe_rms = on_axis_sys.gs.wfe_rms_10e(-9)[0];
            a_wfe_var += wfe_rms * wfe_rms;
            println!("WFE RMS: {}/{}nm", wfe_rms_0, wfe_rms);
            atm.reset()
        }
        a_wfe_rms.push((a_wfe_var / n_sample as f32).sqrt());
    }
    println!("A. WFE RMS: {:?}", a_wfe_rms);
}
