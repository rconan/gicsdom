use bytes::Bytes;
use bytes::BytesMut;
use futures::TryStreamExt;
use log;
use rusoto_core::Region;
use rusoto_lambda::{
    CreateFunctionRequest, DeleteFunctionRequest, FunctionCode, InvocationRequest, Lambda,
    LambdaClient,
};
use rusoto_s3::{GetObjectRequest, ListObjectsV2Request, S3Client, S3};
use serde_pickle as pickle;
use std::boxed::Box;
use std::error::Error;
use std::io::Cursor;
use std::str::FromStr;
use std::sync::Arc;

pub async fn list(
    region: &str,
    bucket: &str,
    prefix: &str,
    suffix: Option<&str>,
) -> Result<Vec<String>, Box<dyn Error>> {
    let region = Region::from_str(region);
    let s3_client = S3Client::new(region.unwrap_or(Region::default()));
    let mut cfd_keys: Vec<String> = vec![];

    let mut request = ListObjectsV2Request {
        bucket: bucket.to_string(),
        prefix: Some(prefix.to_string()),
        ..Default::default()
    };
    loop {
        let response = s3_client.list_objects_v2(request.clone()).await?;
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
                    Ok(obj)
                    //println!("data: {}", obj.len());
                    //data.push(obj);
                }
                Err(e) => {
                    log::error!("Error: {:?}", e);
                    Err(e)
                    //break;
                }
            }
        }))
    }
    let mut data: Vec<Vec<f32>> = vec![];
    let mut data_ok = Ok(vec![]);
    for h in handle {
        match h.await.unwrap() {
            Ok(out) => {
                data.push(out);
            }
            Err(e) => {
                log::error!("Error: {:?}", e);
                data_ok = Err(format!("{}", e));
                break;
            }
        }
    }
    if data_ok.is_ok() {
        data_ok = Ok(data);
    }
    data_ok
}

pub struct AWSLambda {
    client: LambdaClient,
    function_name: String,
}
impl AWSLambda {
    pub fn new(region: &str, function_name: &str) -> Self {
        let aws_region = Region::from_str(region).unwrap_or(Region::default());
        AWSLambda {
            client: LambdaClient::new(aws_region),
            function_name: function_name.to_string(),
        }
    }

    pub async fn create(self, bucket: &str, fun_key: &str, role: &str) -> Result<Self,Box<dyn Error>> {
        let request = CreateFunctionRequest {
            code: FunctionCode {
                s3_bucket: Some(bucket.to_string()),
                s3_key: Some(fun_key.to_string()),
                ..Default::default()
            },
            function_name: self.function_name.clone(),
            handler: "lambda_function.lambda_handler".to_string(),
            memory_size: Some(3000),
            runtime: "python3.6".to_string(),
            role: role.to_string(),
            timeout: Some(180),
            ..Default::default()
        };
        self.client.create_function(request).await?;
        Ok(self)
    }

    pub async fn delete(self) -> Result<Self,Box<dyn Error>> {
        let request = DeleteFunctionRequest {
            function_name: self.function_name.clone(),
            ..Default::default()
        };
        self.client.delete_function(request).await?;
        Ok(self)
    }

    pub async fn invoke(self, payloads: &[String]) -> Result<Self,Box<dyn Error>> {
        let mut request = InvocationRequest {
            function_name: self.function_name.clone(),
            invocation_type: Some("Event".to_string()),
            ..Default::default()
        };
        for payload in payloads {
            request.payload = Some(Bytes::from(payload.to_string()));
            self.client.invoke(request.clone()).await?;
        }
        Ok(self)
    }
}