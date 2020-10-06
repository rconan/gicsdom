use cirrus;
use log::LevelFilter;
use simple_logger::SimpleLogger;
use serde_pickle as pickle;
use std::time::Instant;
use std::fs::File;
use std::io::BufWriter;

#[tokio::main]
async fn main() {
    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .with_module_level("cirrus", LevelFilter::Warn)
        .init()
        .unwrap();

    let cfd_keys = cirrus::list(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_",
        None,
    )
    .await.unwrap();
    let n_cfd_keys = cfd_keys.len();
    println!(
        "CFD keys #: {} ; [{},...,{}]",
        n_cfd_keys,
        cfd_keys[0],
        cfd_keys[n_cfd_keys - 1]
    );

    let n_opd = 10_usize;
    let now = Instant::now();
    let opd = cirrus::load("us-west-2", "gmto.modeling", &cfd_keys[0..n_opd]).await.unwrap();
    println!("Downloaded {} files in {}s",n_opd,now.elapsed().as_secs());

    let file = File::create("OPD.pkl").unwrap();
    let mut writer = BufWriter::with_capacity(1_000_000, file);
    pickle::to_writer(&mut writer, &opd, true).unwrap();
}
