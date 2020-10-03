use bytes::Bytes;
use log;
use rusoto_core::Region;
use rusoto_s3::{ListObjectsV2Request, S3Client, S3};
use serde_json::json;
use std::str::FromStr;

pub async fn list(region: &str, bucket: &str, prefix: &str, suffix: Option<&str>) -> Vec<String> {
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
                    None => break,
                }
            }
            Err(error) => {
                log::error!("Error: {:?}", error);
                break;
            }
        }
    }
    cfd_keys
}
