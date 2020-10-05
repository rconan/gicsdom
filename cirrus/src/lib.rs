use bytes::Bytes;
use bytes::BytesMut;
use futures::TryStreamExt;
use indicatif::{ProgressBar, ProgressStyle};
use log;
use rusoto_core::Region;
use rusoto_lambda::{InvocationRequest, Lambda, LambdaClient};
use rusoto_s3::{GetObjectRequest, ListObjectsV2Request, S3Client, S3};
use serde_json::json;
use serde_pickle as pickle;
use std::io::Cursor;
use std::str::FromStr;
use std::sync::Arc;
use std::thread::sleep;
use std::time::Duration;

pub async fn list(
    region: &str,
    bucket: &str,
    prefix: &str,
    suffix: Option<&str>,
) -> Result<Vec<String>, String> {
    let region = Region::from_str(region);
    let s3_client = S3Client::new(region.unwrap_or(Region::default()));
    let mut cfd_keys: Vec<String> = vec![];

    let mut request = ListObjectsV2Request {
        bucket: bucket.to_string(),
        prefix: Some(prefix.to_string()),
        ..Default::default()
    };
    loop {
        let result = s3_client.list_objects_v2(request.clone()).await;
        match result {
            Ok(response) => {
                let mut objects = response
                    .contents
                    .unwrap()
                    .iter()
                    .map(|x| x.key.as_ref().cloned().unwrap())
                    .filter(|x| suffix.map_or_else(|| true, |v| x.ends_with(v)))
                    .collect::<Vec<String>>();
                log::info!("Keys #: {}", objects.len());
                cfd_keys.append(&mut objects);
                match response.next_continuation_token {
                    Some(token) => request.continuation_token = Some(token),
                    None => break Ok(cfd_keys),
                }
            }
            Err(error) => {
                log::error!("Error: {:?}", error);
                break Err("List failed!".to_string());
            }
        }
    }
}

pub async fn load_serial(region: &str, bucket: &str, keys: &[String]) -> Vec<Vec<f32>> {
    let region = Region::from_str(region);
    let s3_client = S3Client::new(region.unwrap_or(Region::default()));

    let mut data: Vec<Vec<f32>> = vec![];
    for key in keys {
        log::info!("Downloading {}...", key);
        let request = GetObjectRequest {
            bucket: bucket.to_string(),
            key: key.clone(),
            ..Default::default()
        };
        let result = s3_client.get_object(request).await;
        match result {
            Ok(response) => {
                let stream = response.body.unwrap();
                let body = stream
                    .map_ok(|b| BytesMut::from(&b[..]))
                    .try_concat()
                    .await
                    .unwrap();
                let r = Cursor::new(body);
                let obj: Vec<f32> = pickle::from_reader(r).unwrap();
                //println!("data: {}", obj.len());
                data.push(obj);
            }
            Err(error) => {
                log::error!("Error: {:?}", error);
                break;
            }
        }
    }
    data
}
pub async fn load(region: &str, bucket: &str, keys: &[String]) -> Result<Vec<Vec<f32>>, String> {
    let region = Region::from_str(region);
    let s3_client = Arc::new(S3Client::new(region.unwrap_or(Region::default())));

    let mut handle = vec![];
    for key in keys {
        log::info!("Downloading {}...", key);
        let request = GetObjectRequest {
            bucket: bucket.to_string(),
            key: key.clone(),
            ..Default::default()
        };
        let s3_client = s3_client.clone();
        handle.push(tokio::spawn(async move {
            let result = s3_client.get_object(request).await;
            match result {
                Ok(response) => {
                    let stream = response.body.unwrap();
                    let body = stream
                        .map_ok(|b| BytesMut::from(&b[..]))
                        .try_concat()
                        .await
                        .unwrap();
                    let r = Cursor::new(body);
                    let obj: Vec<f32> = pickle::from_reader(r).unwrap();
                    Some(obj)
                    //println!("data: {}", obj.len());
                    //data.push(obj);
                }
                Err(error) => {
                    log::error!("Error: {:?}", error);
                    None
                    //break;
                }
            }
        }))
    }
    let mut data: Vec<Vec<f32>> = vec![];
    let mut data_ok = true;
    for h in handle {
        match h.await.unwrap() {
            Some(out) => {
                data.push(out);
            }
            None => {
                log::warn!("No data!");
                data_ok = false;
                break;
            }
        }
    }
    if data_ok {
        Ok(data)
    } else {
        Err("Corrupted data!".to_string())
    }
}

pub async fn invoke(function_name: &str, keys: Vec<String>) {
    let lambda_client = LambdaClient::new(Region::UsEast2);

    let mut request = InvocationRequest {
        function_name: function_name.to_string(),
        invocation_type: Some("Event".to_string()),
        ..Default::default()
    };

    let pause = Duration::from_secs(60);
    let chunck_size = 250_usize;
    for keys in keys.chunks(chunck_size) {
        //println!("Processing {} of keys",keys.len());
        let bar = ProgressBar::new(chunck_size as u64);
        bar.set_style(
            ProgressStyle::default_bar().template("[{bar:100.cyan/blue} {pos:>5}/{len:5}"),
        );
        for key in keys {
            bar.inc(1);
            let payload = json!({
                "bucket": "gmto.starccm",
                "key": key
            });
            request.payload = Some(Bytes::from(payload.to_string()));
            let result = lambda_client.invoke(request.clone()).await;
            match result {
                Ok(_response) => (), //println!("Response: {:?}", response),
                Err(error) => {
                    println!("Error: {:?}", error);
                    break;
                }
            }
        }
        bar.finish_and_clear();
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(1000);
        pb.set_style(ProgressStyle::default_spinner());
        //println!("Waiting for the 1st wave of lambdas to be executed ...");
        sleep(pause);
        pb.finish_and_clear();
    }
}
