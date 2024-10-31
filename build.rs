use std::{env, path::Path, process::Command};

fn main() {
    let abpos_rel_dir_name = "abPOA-v1.5.3";
    let out_path_str = env::var("OUT_DIR").unwrap();
    let build_out_path = Path::new(&out_path_str);
    // let c_cur_dir = env::current_dir().unwrap();
    // let build_out_path = Path::new(c_cur_dir.to_str().unwrap());

    Command::new("sh")
        .arg("-c")
        .arg(&format!(
            "cp -r {} {}",
            abpos_rel_dir_name,
            build_out_path.to_str().unwrap()
        ))
        .output()
        .expect("cp abpoa_src error");

    let c_src_file_dir = build_out_path.join(abpos_rel_dir_name);
    let c_src_file_dir = &c_src_file_dir;
    let current_dir = env::current_dir().unwrap();

    Command::new("sed")
        .current_dir(c_src_file_dir)
        .arg("-i")
        .arg("s/EXTRA_FLAGS = -Wall -Wno-unused-function -Wno-misleading-indentation -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES# -fno-tree-vectorize/EXTRA_FLAGS = -Wall -Wno-unused-function -Wno-misleading-indentation -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES -fPIC # -fno-tree-vectorize/g")
        .arg("Makefile")
        .output()
        .expect("Failed to modify abpoa makefile.");

    // // OUT_DIR is the build output dir
    // let build_out_path = Path::new(&env::var("OUT_DIR").unwrap());

    Command::new("make")
        .arg("-j8")
        .arg("libabpoa")
        .current_dir(c_src_file_dir.to_str().unwrap())
        .output()
        .expect("build error");

    let obj_filedir = c_src_file_dir.join("lib");
    let obj_filedir = &obj_filedir;

    // Command::new("sh")
    //     .arg("-c")
    //     .arg("rm *.o")
    //     .current_dir(obj_filedir)
    //     .output()
    //     .expect("clean error");

    println!(
        "cargo:rerun-if-changed={}",
        current_dir.to_str().unwrap()
    );
    println!(
        "cargo:rustc-link-search=native={}",
        obj_filedir.to_str().unwrap()
    );
    println!("cargo:rustc-link-lib=static=abpoa");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");
}
