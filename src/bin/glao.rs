use std::time::Instant;
use gicsdom::ceo;

fn ref_pssn() {

    let npx = 512;
    let n_src: usize = 6;
    let mut src = ceo::Source::new(n_src as i32,25.5,npx);
    let zen: Vec<f32> = (0..n_src).map(|x| ceo::Conversion::from_arcmin(x as f32) ).collect();
    let azi = vec![0.0;n_src];
    let mag = vec![0.0;n_src];
    src.build("V",zen,azi,mag);

    let mut gmt = ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(None);
    src.through(&mut gmt).xpupil();
    println!("WFE RMS: {:?}nm",src.wfe_rms_10e(-9)[0]);

    let mut pssn: ceo::GPSSn<ceo::AtmosphereTelescopeError> = ceo::GPSSn::new();
    pssn.build(&mut src);
    println!("PSSn: {}",pssn.reset(&mut src));

    let mut atm = ceo::Atmosphere::new();
    atm.gmt_build(pssn.r0(), pssn.oscale);

    let n_sample = 100;
    let now = Instant::now();
    for _ in 0..n_sample {
        src.through(&mut gmt).xpupil().through(&mut atm);
        pssn.accumulate(&mut src);
        //println!("WFE RMS: {:?}nm",src.wfe_rms_10e(-9));
        //println!("PSSn: {}",pssn.peek(&mut src));
        atm.reset();
    }
    println!("ET: {:.3}mn",now.elapsed().as_millis() as f64*1e-3/60f64);
    println!("PSSn: {}",pssn.reset(&mut src));

}

fn main() {

}
