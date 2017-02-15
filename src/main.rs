#[macro_use]
extern crate serde_derive;

use std::collections::HashMap;
use std::fs::File;
use std::io::{Write, Read};
use std::path::Path;

extern crate clap;
extern crate murmurhash3;
extern crate needletail;
extern crate serde;
extern crate serde_json;

use clap::{App, AppSettings, Arg, SubCommand};
use needletail::fastx::fastx_file;
use needletail::seq::Seq;
use needletail::kmer::normalize;

mod distance;
mod filtering;
mod minhashes;
mod serialization;

use distance::distance;
use filtering::filter_sketch;
use minhashes::MinHashKmers;
use serialization::{JSONSketch, JSONMultiSketch};

macro_rules! add_output_options {
    ($cmd:ident) => {
        $cmd = $cmd.arg(Arg::with_name("output_file")
             .short("o")
             .long("output")
             .help("Output to this file")
             .conflicts_with("std_out")
             .takes_value(true))
        .arg(Arg::with_name("std_out")
             .short("O")
             .long("std-out")
             .help("Output to stdout ('print to terminal')")
             .conflicts_with("output_file"));
    };
}

macro_rules! output_to {
    ($object: ident, $matches: ident) => {
        let output = $matches.value_of("output_file");
        match output {
            None => println!("{}", serde_json::to_string(&$object).unwrap()),
            Some(out_filename) => {
                let mut out = File::create(out_filename).unwrap();
                let _ = out.write_all(&serde_json::to_vec(&$object).unwrap());
            },
        }
    }

}

macro_rules! add_kmer_options {
    ($cmd:ident) => {
        $cmd = $cmd.arg(Arg::with_name("n_hashes")
             .short("n")
             .long("n-hashes")
             .help("How many kmers/hashes to store")
             .takes_value(true)
             .default_value("2000"))
        .arg(Arg::with_name("seed")
             .long("seed")
             .help("Seed murmurhash with this value")
             .takes_value(true)
             .default_value("42"))
        .arg(Arg::with_name("kmer_length")
             .short("k")
             .long("kmer-length")
             .help("Length of kmers to use")
             .takes_value(true)
             .default_value("21"))
        .arg(Arg::with_name("filter")
             .short("f")
             .long("filter")
             .help("Filter kmers with an cumulative abundance lower than this out")
             .takes_value(true)
             .default_value("1.0"))
        .arg(Arg::with_name("oversketch")
             .long("oversketch")
             .help("The amount of extra sketching to do before filtering. This is only a safety to allow sketching e.g. high-coverage files with lots of error-generated uniquemers and should not change the final sketch")
             .requires("filter")
             .takes_value(true)
             .default_value("100"))
        .arg(Arg::with_name("no_strict")
             .short("N")
             .long("no-strict")
             .help("Allow sketching files with fewer kmers than `n_hashes`"));
    }
}


fn main() {
    let mut sketch_command = SubCommand::with_name("sketch")
        .about("Sketch FASTA/Q file(s) into MASH sketches")
        .arg(Arg::with_name("INPUT")
             .help("The file(s) to sketch")
             .multiple(true)
             .required(true));
    add_output_options!(sketch_command);
    add_kmer_options!(sketch_command);

    let mut dist_command = SubCommand::with_name("dist")
        .about("Compute distances between MASH sketches")
        .arg(Arg::with_name("INPUT")
             .help("Sketchfile to make comparisons for")
             // .multiple(true)  // TODO
             .required(true))
        .arg(Arg::with_name("pairwise")
             .short("p")
             .long("pairwise")
             .conflicts_with("queries")
             .help("Calculate distances between all sketches"))
        .arg(Arg::with_name("queries")
             .help("All distances are from these sketches")
             .multiple(true)
             .conflicts_with("pairwise")
             .takes_value(true))
        .arg(Arg::with_name("mash_mode")
             .short("m")
             .long("mash")
             .help("Calculate distances using the same algorithms as Mash"));
    add_output_options!(dist_command);
    add_kmer_options!(dist_command);

    let mut join_command = SubCommand::with_name("join")
        .about("Merge multiple sketch files into a multisketch file")
        .arg(Arg::with_name("error_filter")
             .short("f")
             .long("filter")
             .help("Filter kmers with an abundance lower than this out")
             .takes_value(true)
             .default_value("1.0"));
    add_output_options!(join_command);

    let mut hist_command = SubCommand::with_name("hist")
        .arg(Arg::with_name("INPUT")
             .help("Generate histograms from these file(s)")
             .multiple(true)
             .required(true));
    add_output_options!(hist_command);
    add_kmer_options!(hist_command);

    let mut info_command = SubCommand::with_name("info")
        .arg(Arg::with_name("INPUT")
             .help("Return stats on these file(s)")
             .multiple(true)
             .required(true));
    add_output_options!(info_command);
    add_kmer_options!(info_command);

    let matches = App::new("finch").version("0.1.0")
        .author("Roderick Bovee <roderick@onecodex.com>")
        .about("Work with MASH sketches")
        .setting(AppSettings::VersionlessSubcommands)
        .subcommand(sketch_command)
        .subcommand(dist_command)
        .subcommand(join_command)
        .subcommand(hist_command)
        .subcommand(info_command)
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("sketch") {
        let final_sketch_size = matches.value_of("n_hashes").unwrap().parse::<usize>().expect("n_hashes must be an integer");
        let seed = matches.value_of("seed").unwrap().parse::<u64>().expect("seed must be an integer");
        let kmer_length = matches.value_of("kmer_length").unwrap().parse::<u8>().expect("kmer_length must be an integer < 256");
        let no_strict = matches.is_present("no_strict");

        let (sketch_size, filter) = match matches.value_of("filter") {
            Some(filter) => {
                let oversketch = matches.value_of("oversketch").unwrap().parse::<usize>().unwrap();
                (oversketch * final_sketch_size, filter.parse::<f32>().expect("filter must be a number between 0 and 100") / 100f32)
            },
            None => (final_sketch_size, 0f32),
        };

        let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
        if matches.is_present("output_file") || matches.is_present("std_out") {
            let sketches = mash_files(filenames, sketch_size, final_sketch_size, kmer_length, filter, no_strict, seed);
            match sketches {
                Ok(s) => {
                    output_to!(s, matches);
                },
                Err(e) => panic!(e),
            }
        } else {
            // "sketch in place"
            for filename in &filenames {
                let out_filename = filename.to_string() + ".sk";
                let sketches = mash_files(vec![filename], sketch_size, final_sketch_size, kmer_length, filter, no_strict, seed);
                match sketches {
                    Ok(s) => {
                        let mut out = File::create(out_filename).unwrap();
                        let _ = out.write_all(&serde_json::to_vec(&s).unwrap());
                    },
                    Err(e) => panic!(e),
                }
            }
        }
    } else if let Some(matches) = matches.subcommand_matches("dist") {
        // TODO: handle the multiple files case
        let filename = matches.value_of("INPUT").unwrap();
        let mash_mode = matches.is_present("mash_mode");
        let file = match File::open(filename) {
            Ok(v) => v,
            Err(e) => panic!("Error opening file: {}", e),
        };
        let json: JSONMultiSketch = match serde_json::from_reader(file) {
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


fn mash_files(filenames: Vec<&str>, n_hashes: usize, final_size: usize, kmer_length: u8, filter: f32, no_strict: bool, seed: u64) -> Result<JSONMultiSketch, String> {
    let mut sketches = Vec::with_capacity(filenames.len());
    for filename in &filenames {
        let sketch = mash_file(Path::new(filename), n_hashes, final_size, kmer_length, filter, seed);
        if !no_strict && sketch.len() < final_size {
            return Err(format!("{} had too few kmers ({}) to sketch", filename, sketch.len()));
        }
        sketches.push(sketch);
    }
    Ok(JSONMultiSketch {
        kmer: kmer_length,
        alphabet: String::from("ACGT"),
        preserveCase: false,
        canonical: true,
        sketchSize: final_size as u32,
        hashType: String::from("MurmurHash3_x64_128"),
        hashBits: 64u16,
        hashSeed: 42u64,
        sketches: sketches,
    })
}


fn mash_file(filename: &Path, n_hashes: usize, final_size: usize, kmer_length: u8, filter: f32, seed: u64) -> JSONSketch {
    let mut minhash = MinHashKmers::new(n_hashes, seed);
    let mut seq_len = 0u64;
    fastx_file(filename.to_str().unwrap(), |seq| {
        let norm_seq = normalize(&seq.1, false);
        let mut norm_seq = Seq::new(&norm_seq);
        for kmer in norm_seq.canonical_kmers(kmer_length) {
            minhash.push(kmer);
        }
        seq_len += seq.1.len() as u64;
    }).unwrap();

    let mut hashes = minhash.into_vec();
    let mut filter_stats: HashMap<String, String> = HashMap::new();
    if filter > 0f32 {
        let (filtered_hashes, cutoff, _) = filter_sketch(&hashes, filter, final_size);
        hashes = filtered_hashes;
        filter_stats.insert(String::from("minCopies"), cutoff.to_string());
    }

    // directory should be clipped from filename
    let basename = filename.file_name().unwrap();
    JSONSketch::new(basename.to_str().unwrap(), seq_len, hashes, &filter_stats)
}
