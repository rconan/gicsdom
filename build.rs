extern crate bindgen;

//use std::env;
//use std::path::PathBuf;

fn main() {
    let bindings = bindgen::Builder::default()
        .header("wrapper.hpp")
        .clang_arg("-v")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .whitelist_type("source")
        .whitelist_type("imaging")
        .whitelist_type("shackHartmann")
        .whitelist_type("geometricShackHartmann")
        .whitelist_type("gmt_m1")
        .whitelist_type("gmt_m2")
        .whitelist_type("atmosphere")
        .generate()
        .expect("Unable to generate bindings");
//    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file("src/bindings.rs")
        .expect("Couldn't write bindings!");
    println!("cargo:rustc-link-search=native=CEO/lib/");
    println!("cargo:rustc-link-lib=static=ceo");
    println!("cargo:rustc-link-lib=curl");
    println!("cargo:rustc-link-lib=jsmn");
    println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    println!("cargo:rustc-link-lib=cudart");
    println!("cargo:rustc-link-lib=cublas");
    println!("cargo:rustc-link-lib=cufft");
    println!("cargo:rustc-link-lib=cusparse");
    println!("cargo:rustc-link-lib=curand");
    println!("cargo:include=CEO/include");
    println!("cargo:lib=CEO/lib");
    println!("cargo:rerun-if-changed=wrapper.hpp");
}
