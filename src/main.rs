#[macro_use]
extern crate serde_derive;

use std::fs::File;
use std::io::{Write, Read};
use std::path::Path;

extern crate clap;
extern crate murmurhash3;
extern crate needletail;
extern crate serde;
extern crate serde_json;

use clap::{App, Arg, SubCommand};
use needletail::fastx::fastx_file;
use needletail::kmer::{canonical, has_no_n, normalize};

mod distance;
mod minhashes;
mod serialization;

use distance::distance;
use minhashes::MinHashKmers;
use serialization::JSONSketch;


fn main() {
    let sketch_command = SubCommand::with_name("sketch")
        .about("Sketch FASTA/Q file(s) into MASH sketches")
        .arg(Arg::with_name("INPUT")
             .help("The file(s) to sketch")
             .multiple(true)
             .required(true))
        .arg(Arg::with_name("n_hashes")
             .short("n")
             .long("n-hashes")
             .help("How many kmers/hashes to store")
             .takes_value(true)
             .default_value("2000"))
        .arg(Arg::with_name("kmer_length")
             .short("k")
             .long("kmer-length")
             .help("Length of kmers to use")
             .takes_value(true)
             .default_value("21"))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .help("File to output sketch to")
             .takes_value(true));

    let dist_command = SubCommand::with_name("dist")
        .about("Compute distanes between MASH sketches")
        .arg(Arg::with_name("INPUT")
             .help("Sketchfile to make comparisons for")
             .required(true))
        .arg(Arg::with_name("mash_mode")
             .short("m")
             .long("mash")
             .help("Calculate distances as Mash"));

    let matches = App::new("finch").version("0.1.0")
        .author("Roderick Bovee <roderick@onecodex.com>")
        .about("Work with MASH sketches")
        .subcommand(sketch_command)
        .subcommand(dist_command)
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("sketch") {
        let kmer_length = matches.value_of("kmer_length").unwrap().parse::<u8>().unwrap();
        let n_hashes = matches.value_of("n_hashes").unwrap().parse::<usize>().unwrap();
        let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
        let output = matches.value_of("output");

        let mut sketches = Vec::with_capacity(filenames.len());
        for filename in &filenames {
            sketches.push(mash_file(Path::new(filename), n_hashes, kmer_length));
        }
        let output_data = JSONOutput {
            kmer: kmer_length,
            alphabet: String::from("ACGT"),
            preserveCase: false,
            canonical: true,
            sketchSize: n_hashes as u32,
            hashType: String::from("MurmurHash3_x64_128"),
            hashBits: 64u16,
            hashSeed: 42u16,
            sketches: sketches,
        };
        match output {
            None => println!("{}", serde_json::to_string(&output_data).unwrap()),
            Some(out_filename) => {
                let mut out = File::create(out_filename).unwrap();
                out.write_all(&serde_json::to_vec(&output_data).unwrap());
            },
        };
    } else if let Some(matches) = matches.subcommand_matches("dist") {
        let filename = matches.value_of("INPUT").unwrap();
        let mash_mode = matches.is_present("mash_mode");
        let file = match File::open(filename) {
            Ok(v) => v,
            Err(e) => panic!("Error opening file: {}", e),
        };
        let json: JSONOutput = match serde_json::from_reader(file) {
            Ok(v) => v,
            Err(e) => panic!("Error parsing file: {}", e),
        };
        for raw_sketch1 in json.sketches.iter() {
            let sketch1 = &raw_sketch1.get_kmers().unwrap();
            for raw_sketch2 in json.sketches.iter() {
                if raw_sketch1 == raw_sketch2 {
                    continue;
                }
                let sketch2 = &raw_sketch2.get_kmers().unwrap();
                let distance_metrics = distance(&sketch1, &sketch2, mash_mode).unwrap();
                let distance = distance_metrics.0;
                let jaccard = distance_metrics.1;
                let common = distance_metrics.2;
                let total = distance_metrics.3;
                println!(
                    "Distance from {} to {}: {} (J{} {}/{})",
                    raw_sketch1.name, raw_sketch2.name, distance, jaccard, common, total
                );
            }
        }
    }
}


#[derive(Serialize, Deserialize)]
struct JSONOutput {
    kmer: u8,
    alphabet: String,
    preserveCase: bool,
    canonical: bool,
    sketchSize: u32,
    hashType: String,
    hashBits: u16,
    hashSeed: u16,
    sketches: Vec<JSONSketch>,
}


fn mash_file(filename: &Path, n_hashes: usize, kmer_length: u8) -> JSONSketch {
    let mut minhash = MinHashKmers::new(n_hashes);
    let mut seq_len = 0u64;
    fastx_file(filename.to_str().unwrap(), |seq| {
        let norm_seq = normalize(&seq.1, false);
        if !has_no_n(&norm_seq) {
            // TODO: write a "no N" kmer iterator (with normalization)
            for kmer in norm_seq.windows(kmer_length as usize) {
                if !has_no_n(kmer) {
                    continue;
                }
                minhash.push(&canonical(kmer));
            }
        } else {
            for kmer in norm_seq.windows(kmer_length as usize) {
                minhash.push(&canonical(kmer));
            }
        }
        seq_len += norm_seq.len() as u64;
    }).unwrap();

    let hashes = minhash.into_vec();
    // directory should be clipped from filename
    let basename = filename.file_name().unwrap();
    JSONSketch::new(basename.to_str().unwrap(), seq_len, hashes)
}
