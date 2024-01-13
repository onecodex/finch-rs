use std::io::{BufReader, Cursor};
use std::process::Command;

use assert_cmd::prelude::*;
use predicates::prelude::predicate;

use finch::serialization::{read_finch_file, read_mash_file, Sketch};

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
fn finch_sketch_bin() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("finch")?;
    cmd.arg("sketch")
        .args(&["--n-hashes", "10"])
        .arg("-b")
        .arg("-O")
        .arg("tests/data/query.fa");
    cmd.assert().success();

    let output = Cursor::new(cmd.output().unwrap().stdout);
    let mut buf_reader = BufReader::new(output);
    let sketch: Vec<Sketch> = read_finch_file(&mut buf_reader)?;
    assert_eq!(sketch.len(), 1);
    assert_eq!(sketch[0].sketch_params.k(), 21);
    assert_eq!(sketch[0].sketch_params.expected_size(), 10);
    assert_eq!(sketch[0].hashes.len(), 10);
    Ok(())
}

#[test]
fn finch_sketch_msh() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("finch")?;
    cmd.arg("sketch")
        .args(&["--n-hashes", "10"])
        .arg("-B")
        .arg("-O")
        .arg("tests/data/query.fa");
    cmd.assert().success();

    let output = Cursor::new(cmd.output().unwrap().stdout);
    let mut buf_reader = BufReader::new(output);
    let sketch: Vec<Sketch> = read_mash_file(&mut buf_reader)?;
    assert_eq!(sketch.len(), 1);
    assert_eq!(sketch[0].sketch_params.k(), 21);
    // mash doesn't set an expected size?
    // assert_eq!(sketch[0].sketch_params.expected_size(), 10);
    assert_eq!(sketch[0].hashes.len(), 10);
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

    let kmers = &sketch["sketches"][0]["kmers"];

    assert_eq!(kmers[0], "ATGCTAGCTACGTAACGTCGC");
    assert_eq!(kmers[1], "CAGTCGATCGATCGTAGCTGA");
    assert_eq!(kmers[2], "CTCAGATGCTGAGCCGGTCTA");
    assert_eq!(kmers[3], "GCTAGCTAGCATCGCTAGCTA");
    assert_eq!(kmers[4], "GACTAGCTAGCTAGCTAGCGA");
    assert_eq!(kmers[5], "CGCTAGCTACGATCGATCGAC");
    assert_eq!(kmers[6], "TAATTTATACGGGCCTATTAA");
    assert_eq!(kmers[7], "GCATCAGCTAGCATCGCTGTA");
    assert_eq!(kmers[8], "AGCCGGTCTACTACTACACAT");
    assert_eq!(kmers[9], "AAGGCCTAACTTAATAGGCCC");

    //assert_eq!(sketch["max_hash"], usize::max_value() / 1000);
    assert_eq!(sketch["hashSeed"], 0);

    Ok(())
}
