#[macro_use]
extern crate serde_derive;

use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::File;
use std::io::{Write, Read};
use std::path::Path;

extern crate clap;
extern crate murmurhash3;
extern crate needletail;
extern crate serde;
extern crate serde_json;

use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use needletail::fastx::fastx_cli;
use needletail::seq::Seq;
use needletail::kmer::normalize;

mod distance;
mod filtering;
mod minhashes;
mod serialization;

use distance::distance;
use filtering::{filter_sketch, filter_strands, hist};
use minhashes::MinHashKmers;
use serialization::{JSONSketch, JSONMultiSketch};

const FINCH_EXT: &'static str = ".sk";

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
             .default_value("0"))
        .arg(Arg::with_name("kmer_length")
             .short("k")
             .long("kmer-length")
             .help("Length of kmers to use")
             .takes_value(true)
             .default_value("21"))
        .arg(Arg::with_name("filter")
             .short("f")
             .long("filter")
             .help("Filter out kmers with an cumulative abundance lower than this (error filtering)")
             .takes_value(true)
             .default_value("1.0"))
        .arg(Arg::with_name("strand_filter")
             .long("strand-filter")
             .help("Filter out kmers with a canonical kmer percentage lower than this (adapter filtering)")
             .takes_value(true)
             .default_value("0.1"))
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
             .help("Sketchfile(s) to make comparisons for")
             .multiple(true)
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
        .subcommand(hist_command)
        .subcommand(info_command)
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("sketch") {
        let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
        if matches.is_present("output_file") || matches.is_present("std_out") {
            let sketches = open_mash_files(filenames, matches);
            match sketches {
                Ok(s) => {
                    output_to!(s, matches);
                },
                Err(e) => panic!(e),
            }
        } else {
            // "sketch in place"
            for filename in &filenames {
                let out_filename = filename.to_string() + FINCH_EXT;
                let sketches = open_mash_files(vec![filename], matches);
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
        let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
        let mash_mode = matches.is_present("mash_mode");

        let all_sketch_objs = open_mash_files(filenames, matches).unwrap();
        let all_sketches = all_sketch_objs.sketches;

        let mut query_pos: HashSet<usize> = HashSet::with_capacity(all_sketches.len());
        if matches.is_present("pairwise") {
            query_pos = (0..all_sketches.len()).collect();
        } else if matches.is_present("queries") {
            let query_names: HashSet<String> = matches.values_of("queries").unwrap().map(|s| s.to_string()).collect();

            for (i, sketch) in all_sketches.iter().enumerate() {
                if query_names.contains(&sketch.name) {
                    query_pos.insert(i);
                }
            }
        } else {
            query_pos.insert(0);
        }
        let mut distances = Vec::new();
        for raw_sketch1 in all_sketches.iter()
                                       .enumerate()
                                       .filter(|&(i, _)| query_pos.contains(&i))
                                       .map(|(_, s)| s)
        {
            let sketch1 = &raw_sketch1.get_kmers().unwrap();
            for raw_sketch2 in all_sketches.iter() {
                if raw_sketch1 == raw_sketch2 {
                    continue;
                }
                let sketch2 = &raw_sketch2.get_kmers().unwrap();
                let distance = distance(&sketch1, &sketch2, &raw_sketch1.name, &raw_sketch2.name, mash_mode).unwrap();
                distances.push(distance);
            }
        }
        output_to!(distances, matches);
    } else if let Some(matches) = matches.subcommand_matches("hist") {
        let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();

        let all_sketch_objs = open_mash_files(filenames, matches).unwrap();
        let all_sketches = all_sketch_objs.sketches;

        let mut hist_map: HashMap<&str, Vec<u64>> = HashMap::new();
        for sketch in all_sketches.iter() {
            hist_map.insert(&sketch.name, hist(&sketch.get_kmers().unwrap()));
        }
        output_to!(hist_map, matches);
    }
}


fn mash_files(filenames: Vec<&str>, n_hashes: usize, final_size: usize, kmer_length: u8, filters: &mut FilterParams, no_strict: bool, seed: u64) -> Result<JSONMultiSketch, String> {
    let mut sketches = Vec::with_capacity(filenames.len());
    for filename in &filenames {
        let mut minhash = MinHashKmers::new(n_hashes, seed);
        let mut seq_len = 0u64;
        let path = Path::new(filename);
        fastx_cli(path.to_str().ok_or("Couldn't make path into string")?, |seq_type| {
            // disable filtering for FASTA files unless it was explicitly specified
            if seq_type == "FASTA" {
                if !filters.abun_filter_specified {
                    filters.abun_filter = 0f32;
                }
                if !filters.strand_filter_specified {
                    filters.strand_filter = 0f32;
                }
            }
        }, |seq| {
            let norm_seq = normalize(&seq.1, false);
            let mut norm_seq = Seq::new(&norm_seq);
            for (kmer, is_rev_complement) in norm_seq.canonical_kmers(kmer_length) {
                let rc_count = match is_rev_complement {
                    true => 1u8,
                    false => 0u8,
                };
                minhash.push(kmer, rc_count);
            }
            seq_len += seq.1.len() as u64;
        }).map_err(|e| e.to_string())?;

        let mut hashes = minhash.into_vec();
        let mut filter_stats: HashMap<String, String> = HashMap::new();
        if filters.abun_filter > 0f32 {
            let (filtered_hashes, cutoff) = filter_sketch(&hashes, filters.abun_filter);
            hashes = filtered_hashes;
            filter_stats.insert(String::from("abunFilter"), filters.abun_filter.to_string());
            filter_stats.insert(String::from("minCopies"), cutoff.to_string());
        }

        if filters.strand_filter > 0f32 {
            hashes = filter_strands(&hashes, filters.abun_filter);
            filter_stats.insert(String::from("strandFilter"), filters.strand_filter.to_string());
        }
        hashes.truncate(final_size);

        // directory should be clipped from filename
        let basename = path.file_name().ok_or("Couldn't get filename from path")?;
        let sketch = JSONSketch::new(basename.to_str().ok_or("Couldn't make filename into string")?, seq_len, hashes, &filter_stats);

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
        hashSeed: seed,
        sketches: sketches,
    })
}


fn open_mash_files(filenames: Vec<&str>, matches: &ArgMatches) -> Result<JSONMultiSketch, String> {
    let mut filename_iter = filenames.iter();
    let filename = filename_iter.next().ok_or("At least one filename must be specified")?;

    let mut all_sketch_objs = open_mash_file(filename, matches, None)?;
    for filename in filename_iter {
        let new_sketch = open_mash_file(filename, matches, Some(&all_sketch_objs))?;
        all_sketch_objs.sketches.extend_from_slice(&new_sketch.sketches);
    }
    Ok(all_sketch_objs)
}

struct FilterParams {
    abun_filter: f32,
    abun_filter_specified: bool,
    strand_filter: f32,
    strand_filter_specified: bool,
}

fn open_mash_file(filename: &str, matches: &ArgMatches, default_sketch: Option<&JSONMultiSketch>) -> Result<JSONMultiSketch, String> {
    let final_sketch_size = matches.value_of("n_hashes").unwrap().parse::<usize>().map_err(|_| "n_hashes must be an integer")?;
    let seed = matches.value_of("seed").unwrap().parse::<u64>().map_err(|_| "seed must be an integer")?;
    let kmer_length = matches.value_of("kmer_length").unwrap().parse::<u8>().map_err(|_| "kmer_length must be an integer < 256")?;

    let no_strict = matches.is_present("no_strict");

    let raw_abun_filter = matches.value_of("filter").unwrap().parse::<f32>().map_err(|_| format!("filter must be a percent between 0 and {}", 100f32 / kmer_length as f32))?;

    let raw_strand_filter = matches.value_of("strand_filter").unwrap().parse::<f32>().map_err(|_| format!("strand-filter must be ratio between 0 and 1"))?;

    // note: is_present returns true while occurrences_of correctly is 0
    let mut filters = FilterParams {
        abun_filter: kmer_length as f32 * raw_abun_filter / 100f32,
        abun_filter_specified: matches.occurrences_of("filter") > 0,
        strand_filter: raw_strand_filter,
        strand_filter_specified: matches.occurrences_of("strand_filter") > 0,
    };

    let oversketch = matches.value_of("oversketch").unwrap().parse::<usize>().map_err(|_| "bad value for oversketch")?;
    let sketch_size = final_sketch_size * oversketch;

    // if the file isn't a sketch file, we try to sketch it and pass the sketches back
    if !filename.ends_with(FINCH_EXT) {
        return match default_sketch {
            Some(s) => mash_files(vec![filename], sketch_size / final_sketch_size * s.sketchSize as usize, s.sketchSize as usize, s.kmer, &mut filters, no_strict, s.hashSeed),
            None => mash_files(vec![filename], sketch_size, final_sketch_size, kmer_length, &mut filters, no_strict, seed),
        };
    }

    // otherwise we just open the file and return the sketches
    let file = File::open(filename).map_err(|_|
        format!("Error opening {}", &filename)
    )?;
    let json: JSONMultiSketch = serde_json::from_reader(file).map_err(|_|
        format!("Error parsing {}", &filename)
    )?;

    // sanity checking to make sure different files have comparable hashing parameters
    if let Some(s) = default_sketch {
        // kmer, hashType, hashSeed, and hashBits must be same
        if s.kmer != json.kmer {
            return Err(format!("{} has a different kmer length ({}) from others ({})", filename, s.kmer, json.kmer));
        } else if s.hashType != json.hashType {
            return Err(format!("{} used a different hash ({}) from others ({})", filename, s.hashType, json.hashType));
        } else if s.hashSeed != json.hashSeed {
            return Err(format!("{} had a different hash seed ({}) from others ({})", filename, s.hashSeed, json.hashSeed));
        } else if s.hashBits != json.hashBits {
            return Err(format!("{} used a different length hash ({}) from others ({})", filename, s.hashBits, json.hashBits));
        }
    }
    Ok(json)
}
