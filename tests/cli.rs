use anyhow::Result;
use assert_cmd::Command;
use pretty_assertions::assert_eq;
use rand::{distributions::Alphanumeric, Rng};
use std::{fs, path::Path};
use tempfile::NamedTempFile;

const GEN_COMP: &str = "generate_comparison_data";
const SEQUENCES_3: &str = "tests/inputs/sequences_3.csv";
const SEQUENCES_1K: &str = "tests/inputs/sequences_2001.csv";
const OUTPUT_1K: &str = "tests/outputs/data_2001.csv";
const OUTPUT_3: &str = "tests/outputs/data_3.csv";

// --------------------------------------------------
fn random_string() -> String {
    rand::thread_rng()
        .sample_iter(&Alphanumeric)
        .take(7)
        .map(char::from)
        .collect()
}

// --------------------------------------------------
fn gen_nonexistent_file() -> String {
    loop {
        let filename = random_string();
        if fs::metadata(&filename).is_err() {
            return filename;
        }
    }
}

// --------------------------------------------------
#[test]
fn gen_comp_fails_no_args() {
    let mut cmd = Command::cargo_bin(GEN_COMP).unwrap();
    cmd.assert().failure();
}

// --------------------------------------------------
#[test]
fn gen_comp_fails_bad_file() {
    let mut cmd = Command::cargo_bin(GEN_COMP).unwrap();
    let bad = gen_nonexistent_file();
    cmd.args(["-s", &bad]).assert().failure();
}

// --------------------------------------------------
fn run_gen_comp(seq_file: &str, expected_file: &str) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_str().unwrap();
    let args = &["-s", seq_file, "-o", &outpath];
    let output = Command::cargo_bin(GEN_COMP)?.args(args).output().unwrap();

    // The command succeeds
    assert!(output.status.success());

    // The expected output file exists
    assert!(Path::new(&outpath).exists());

    // The contents of the output file match the expected value
    let expected = fs::read_to_string(expected_file)?;
    let actual = fs::read_to_string(&outpath)?;
    assert_eq!(actual, expected);

    Ok(())
}

// --------------------------------------------------
#[test]
fn run_gen_comp_3() -> Result<()> {
    run_gen_comp(SEQUENCES_3, OUTPUT_3)
}

// --------------------------------------------------
#[test]
fn run_gen_comp_1k() -> Result<()> {
    run_gen_comp(SEQUENCES_1K, OUTPUT_1K)
}
