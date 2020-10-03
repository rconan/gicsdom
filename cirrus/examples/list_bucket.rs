use cirrus;
use log::LevelFilter;
use simple_logger::SimpleLogger;
use std::env;

// Usage:
// cargo run --example list_bucket -- us-east-2  gmto.starccm csv.gz
// cargo run --example list_bucket -- us-west-2  gmto.modeling 

#[tokio::main]
async fn main() {

    SimpleLogger::new()
        .with_level(LevelFilter::Off)
        .with_module_level("cirrus", LevelFilter::Info)
        .init()
        .unwrap();

    let mut args = env::args().skip(1);
    let region = args.next().unwrap();
    let bucket = args.next().unwrap();
    let suffix = args.next().unwrap_or("none".to_owned());
    let q = if suffix=="none" { None } else { Some(suffix.as_str()) };

    let cfd_keys = cirrus::list(
        &region,
        &bucket,
        "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_",
        q,
    )
    .await;
    let n_cfd_keys = cfd_keys.len();
    println!(
        "CFD keys #: {} ; [{},...,{}]",
        n_cfd_keys,
        cfd_keys[0],
        cfd_keys[n_cfd_keys - 1]
    );

    /*
    let cfd_keys = list(
        Region::UsWest2,
        "gmto.modeling",
        "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_",
        None,
    )
    .await;
    let n_cfd_keys = cfd_keys.len();
    println!(
        "CFD keys #: {} ; [{},...,{}]",
        n_cfd_keys,
        cfd_keys[0],
        cfd_keys[n_cfd_keys - 1]
    );
     */
}
