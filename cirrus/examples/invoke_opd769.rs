use cirrus;
use serde::{Deserialize, Serialize};
use serde_json::json;
use serde_yaml;
use std::fs::File;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CfdCases {
    #[serde(rename = "baseline 2020")]
    baseline_2020: Vec<String>,
}

#[tokio::main]
async fn main() {
    let file =
        File::open("/home/rconan/Dropbox/Documents/GMT/CFD/Python/CFD/CFD_CASES.yaml").unwrap();
    let cfd_cases: CfdCases = serde_yaml::from_reader(file).unwrap();
    println!("CFD CASES: {:?}", cfd_cases);

    let lambda_name = "CFD2OPD769";
    let mut handle = vec![];

    for cfd_case in cfd_cases.baseline_2020.iter().take(60) {
        let folder = format!("Baseline2020/{}/OPDData_OPD_Data_", cfd_case);
        let cfd_keys = cirrus::list("us-east-2", "gmto.starccm", &folder, Some(".csv.gz"))
            .await
            .unwrap();
        /*
        let n_cfd_keys = cfd_keys.len();
        println!(
            "CFD keys #: {} ; [{},...,{}]",
            n_cfd_keys,
            cfd_keys[0],
            cfd_keys[n_cfd_keys - 1]
        );
         */
        let payloads = cfd_keys
            .iter()
            .map(|x| {
                let p = json!({
                    "bucket": "gmto.starccm",
                    "key": x
                });
                p.to_string()
            })
            .collect::<Vec<String>>();
        for (k, ppayloads) in payloads.chunks(1000).enumerate() {
            let t_payloads = ppayloads.to_vec();
            let t_cfd_case = String::from(cfd_case);
            handle.push(tokio::spawn(async move {
                let function_name = format!("{}_{}_{}", lambda_name, t_cfd_case, k);
                cirrus::AWSLambda::new("us-east-2", &function_name).delete().await;
                    /*.create(
                        "gmto.starccm",
                        "CFD769.zip",
                        "arn:aws:iam::378722409401:role/EC2LambdaRole",
                    )
                    .await
                    .invoke(t_payloads.as_slice())
                    .await;*/
            }));
        }
    }
    for h in handle {
        h.await.unwrap();
    }
    /*
    let cfd_keys = cirrus::list(
        "us-west-2",
        "gmto.modeling",
        "Baseline2020/b2019_0z_0az_cd_12ms/OPDData_OPD_Data_",
        None,
    )
    .await
    .unwrap();
    let n_cfd_keys = cfd_keys.len();
    println!(
        "CFD keys #: {} ; [{},...,{}]",
        n_cfd_keys,
        cfd_keys[0],
        cfd_keys[n_cfd_keys - 1]
    );
    */
}
