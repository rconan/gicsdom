use std::f64;

use gicsdom::ceo;
use gicsdom::DomeSeeing;

fn std(opd: Vec<f64>) -> f64 {
    let iter = opd.iter().filter(|x| !x.is_nan());
    let n = iter.clone().fold(0, |s, x| s + 1) as f64;
    let m = iter.clone().fold(0.0, |s, x| s + x) / n;
    let s2 = iter.clone().fold(0.0, |s, x| s + (x - m).powi(2)) / n;
    1e9 * s2.sqrt()
}

fn main() {
    let mut domeseeing = DomeSeeing::new(0, 0, "cd", 12, Some(5.0));
    domeseeing.list();
    let opd = domeseeing.load_at(0.0);
    let iter = opd.iter().filter(|x| !x.is_nan());
    let n = iter.clone().fold(0, |s, x| s + 1) as f64;
    let m = iter.clone().fold(0.0, |s, x| s + x) / n;
    let s2 = iter.clone().fold(0.0, |s, x| s + (x - m).powi(2)) / n;
    println!(
        "l={} ; n={} ; m={} nm ; s={} micron",
        opd.len(),
        n,
        m * 1e9,
        1e6 * s2.sqrt()
    );

    let mut gmt = ceo::Gmt::new();
    gmt.build(0, None);
    let mut src = ceo::Source::new(1, 25.5, 769);
    src.build("V", vec![0.0], vec![0.0], vec![0.0]);
    src.through(&mut gmt).xpupil();
    println!("WFE RMS: {:?}nm", src.wfe_rms_10e(-9));

    let mut gopd = ceo::CuFloat::new();
    gopd.malloc(769 * 769);

    for _ in 0..40 {
        let opd = domeseeing.next().unwrap();
        let mut opd_f32: Vec<f32> = opd
            .iter()
            .map(|&x| if x.is_nan() { 0f32 } else { x as f32 })
            .collect();
        print!("WFE RMS: {:8.3}nm", std(opd));
        gopd.up(&mut opd_f32);
        println!(
            " ; {:8.3}nm",
            src.through(&mut gmt)
                .xpupil()
                .add(&mut gopd)
                .wfe_rms_10e(-9)[0]
        );
    }

    /*
    let mut npz = NpzReader::new(
        File::open("/home/rconan/Downloads/OPDData_OPD_Data_5.000000e+02.npz").unwrap(),
    )
    .unwrap();
    let names = npz.names().unwrap();
    println!("{:?}", names);
    let opd: Array1<f64> = npz.by_name("opd.npy").unwrap();
    println!("{:?}", opd.shape());
    let iter = opd.iter().filter(|x| !x.is_nan());
    let n = iter.clone().fold(0, |s, x| s + 1) as f64;
    let m = iter.clone().fold(0.0, |s, x| s + x) / n;
    let s2 = iter.clone().fold(0.0, |s, x| s + (x - m).powi(2)) / n;
    println!("n={} ; m={} ; s={}", n, m * 1e6, 1e6 * s2.sqrt());

    let s3 = S3::new(Region::UsEast2);
    let bucket = Present(Bucket::from_name("gmto.starccm".to_owned()));
    let object = Object::from_key(
        &bucket,
        "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_5.000000e+02.npz".to_owned(),
    );
    let mut body = Vec::new();

    let now = Instant::now();
    s3.get_body(&Present(object))
        .expect("object body")
        .read_to_end(&mut body)
        .unwrap();
    let r = Cursor::new(body);
    let mut npz = NpzReader::new(r).unwrap();
    println!("Loading in {}s", now.elapsed().as_secs());

    let names = npz.names().unwrap();
    println!("{:?}", names);
    let opd: Array1<f64> = npz.by_name("opd.npy").unwrap();
    println!("{:?}", opd.shape());
    let iter = opd.iter().filter(|x| !x.is_nan());
    let n = iter.clone().fold(0, |s, x| s + 1) as f64;
    let m = iter.clone().fold(0.0, |s, x| s + x) / n;
    let s2 = iter.clone().fold(0.0, |s, x| s + (x - m).powi(2)) / n;
    println!("n={} ; m={} ; s={}", n, m * 1e6, 1e6 * s2.sqrt());

    let o = s3.list_objects(
        &bucket,
        "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_".to_owned(),
    );
    let keys: Vec<String> = o
        .map(|x| x.unwrap().key().to_owned())
        .filter(|x| x.ends_with(".npz"))
        .collect();
    let n = keys.len();
    println!("# of keys: {}", n);
    println!("# of keys[0]: {}", keys[0].trim_end_matches(".npz"));
    let f = |s: &str| -> f64 {
        s.trim_end_matches(".npz")
            .split("/")
            .last()
            .unwrap()
            .split("_")
            .last()
            .unwrap()
            .parse::<f64>()
            .unwrap()
    };

    let opd_time: Vec<f64> = keys.iter().map(|x| f(&x) - f(&keys[0])).collect();
    println!("Time range: [{};{}]", opd_time[0], opd_time[n - 1]);

    let t = 120.01;
    let q = opd_time.iter().map(|x| (x-t).powi(2)).fold(f64::INFINITY,|a,x| f64::min(a,x));
    println!("{}",q);
    println!("{:?}",opd_time.iter().position(|&x|  (x-t).powi(2)==q));
    */
}
