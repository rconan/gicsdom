use ceo;

fn main() {
    let mut gmt= ceo::Gmt::new();
    gmt.build_m1("bending modes", 0).build_m2(Some(0));
    let mut gs = ceo::Source::new(1, 25.5, 512);
    gs.build("Vs", vec![0f32], vec![0f32], vec![0f32] );
    gs.through(&mut gmt).xpupil();
    let mut pssn: ceo::PSSn<ceo::pssn::TelescopeError> = ceo::PSSn::new();
    pssn.build(&mut gs);
    pssn.integrate(&mut  gs);
    println!("WFE RMS: {}; PSSN: {}",gs.wfe_rms_10e(-9)[0], pssn.peek());
}
