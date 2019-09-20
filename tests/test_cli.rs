use std::io::{BufReader, Cursor};
use std::process::Command;

use assert_cmd::prelude::*;
use predicates::prelude::predicate;

use finch::serialization::{read_mash_file, MultiSketch};

#[test]
fn file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("finch")?;
    cmd.arg("sketch").arg("test/file/doesnt/exist");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn finch_sketch() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("finch")?;
    cmd.arg("sketch")
        .args(&["--n-hashes", "10"])
        .args(&["-O"])
        .arg("tests/data/query.fa");
    cmd.assert().success();

    let output = Cursor::new(cmd.output().unwrap().stdout);
    let sketch: serde_json::Value = serde_json::from_reader(output)?;
    assert_eq!(sketch["kmer"], 21);
    assert_eq!(sketch["alphabet"], "ACGT");
    assert_eq!(sketch["sketchSize"], 10);
    assert_eq!(sketch["hashSeed"], 0);

    Ok(())
}

#[test]
fn finch_sketch_msh() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("finch")?;
    cmd.arg("sketch")
        .args(&["--n-hashes", "10"])
        .arg("-b")
        .arg("-O")
        .arg("tests/data/query.fa");
    cmd.assert().success();

    let output = Cursor::new(cmd.output().unwrap().stdout);
    let mut buf_reader = BufReader::new(output);
    let sketch: MultiSketch = read_mash_file(&mut buf_reader)?;
    assert_eq!(sketch.kmer, 21);
    assert_eq!(sketch.alphabet, "ACGT");
    assert_eq!(sketch.hash_seed, 0);
    //mash doesn't save this info...
    //  assert_eq!(sketch.sketchSize, 10);

    Ok(())
}

#[test]
fn finch_sketch_scaled() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("finch")?;
    cmd.arg("sketch")
        .args(&["--n-hashes", "10"])
        .args(&["--sketch-type", "scaled"])
        .args(&["--scale", ".001"])
        .arg("tests/data/query.fa")
        .arg("-O");
    cmd.assert().success();

    let output = Cursor::new(cmd.output().unwrap().stdout);
    let sketch: serde_json::Value = serde_json::from_reader(output)?;
    assert_eq!(sketch["kmer"], 21);
    assert_eq!(sketch["alphabet"], "ACGT");
    assert_eq!(sketch["sketchSize"], 10);
    //assert_eq!(sketch["max_hash"], usize::max_value() / 1000);
    assert_eq!(sketch["hashSeed"], 0);

    Ok(())
}
