#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate clap;

use std::fs::File;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

extern crate murmurhash3;
extern crate needletail;
extern crate serde;
extern crate serde_json;

use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use needletail::fastx::fastx_cli;

mod distance;
mod filtering;
mod minhashes;
mod serialization;
mod statistics;

use distance::distance;
use filtering::{FilterParams, filter_sketch};
use minhashes::MinHashKmers;
use serialization::{JSONSketch, JSONMultiSketch, SketchDistance};
use statistics::{hist, cardinality};

const FINCH_EXT: &'static str = ".sk";

macro_rules! add_output_options {
    ($cmd:ident) => {
        $cmd = $cmd.arg(Arg::with_name("output_file")
             .short("o")
             .long("output")
             .help("Output to this file")
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
             .default_value("1000"))
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
        .arg(Arg::with_name("no_filter")
             .long("no-filter")
             .conflicts_with("filter")
             .help("Disable filtering (default for FASTA)"))
        .arg(Arg::with_name("filter")
             .short("f")
             .long("filter")
             .help("Enable filtering (default for FASTQ)"))
        .arg(Arg::with_name("min_abun_filter")
             .long("min-abun-filter")
             .help("Kmers must have at least this coverage to be included")
             .takes_value(true))
        .arg(Arg::with_name("max_abun_filter")
             .long("max-abun-filter")
             .help("Kmers must have a coverage under this to be included")
             .takes_value(true))
        .arg(Arg::with_name("strand_filter")
             .long("strand-filter")
             .help("Filter out kmers with a canonical kmer percentage lower than this (adapter filtering)")
             .takes_value(true)
             .default_value("0.1"))
        .arg(Arg::with_name("err_filter")
             .long("err-filter")
             .help("Dynamically determine a minimum coverage threshold for filtering from the kmer count histogram using an assumed error rate percentage")
             .takes_value(true)
             .default_value("1"))
        .arg(Arg::with_name("oversketch")
             .long("oversketch")
             .help("The amount of extra sketching to do before filtering. This is only a safety to allow sketching e.g. high-coverage files with lots of error-generated uniquemers and should not change the final sketch")
             .takes_value(true)
             .default_value("200"))
        .arg(Arg::with_name("no_strict")
             .short("N")
             .long("no-strict")
             .help("Allow sketching files with fewer kmers than `n_hashes`"));
    }
}


fn main() {
    let mut sketch_command = SubCommand::with_name("sketch")
        .about("Create sketches from FASTA/Q file(s)")
        .arg(Arg::with_name("INPUT")
             .help("The file(s) to sketch")
             .multiple(true)
             .required(true));
    add_output_options!(sketch_command);
    add_kmer_options!(sketch_command);

    let mut dist_command = SubCommand::with_name("dist")
        .about("Compute distances between sketches")
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
             .short("q")
             .long("queries")
             .help("All distances are from these sketches (sketches must be in the first file)")
             .multiple(true)
             .conflicts_with("pairwise")
             .takes_value(true))
        .arg(Arg::with_name("max_distance")
             .short("d")
             .long("max-dist")
             .help("Only report distances under this threshold")
             .default_value("1.0")
             .takes_value(true))
        .arg(Arg::with_name("mash_mode")
             .short("m")
             .long("mash")
             .help("Calculate distances using the same algorithms as Mash"));
    add_output_options!(dist_command);
    add_kmer_options!(dist_command);

    let mut hist_command = SubCommand::with_name("hist")
        .about("Display histograms of kmer abundances")
        .arg(Arg::with_name("INPUT")
             .help("Generate histograms from these file(s)")
             .multiple(true)
             .required(true));
    add_output_options!(hist_command);
    add_kmer_options!(hist_command);

    let mut info_command = SubCommand::with_name("info")
        .about("Display basic statistics")
        .arg(Arg::with_name("INPUT")
             .help("Return stats on these file(s)")
             .multiple(true)
             .required(true));
    add_output_options!(info_command);
    add_kmer_options!(info_command);

    let matches = App::new("finch").version(crate_version!())
        .author("Roderick Bovee & One Codex <roderick@onecodex.com>")
        .about("Tool for working with genomic MinHash sketches")
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(sketch_command)
        .subcommand(dist_command)
        .subcommand(hist_command)
        .subcommand(info_command)
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("sketch") {
        if matches.is_present("output_file") || matches.is_present("std_out") {
            let sketches = parse_all_mash_files(matches).unwrap();
            output_to!(sketches, matches);
        } else {
            // "sketch in place"
            parse_mash_files(matches, |multisketch, filename| {
                let out_filename = filename.to_string() + FINCH_EXT;
                let mut out = File::create(out_filename).unwrap();
                let _ = out.write_all(&serde_json::to_vec(&multisketch).unwrap());
            }).unwrap();
        }
    } else if let Some(matches) = matches.subcommand_matches("dist") {
        let mash_mode = matches.is_present("mash_mode");

        let max_dist = matches.value_of("max_distance").unwrap().parse::<f64>().map_err(|_| {
            format!("max-distance must be a number")
        }).and_then(|r| {
            if 0f64 <= r && r <= 1f64 {
                return Ok(r);
            }
            Err(format!("max-distance must be 0 and 1"))
        }).unwrap();

        // some special logic if we're doing a pairwise comparison
        if matches.is_present("pairwise") {
            let mut distances = Vec::new();
            let all_sketches = parse_all_mash_files(matches).unwrap();
            distances.extend(calc_sketch_distances(&all_sketches.sketches, &all_sketches.sketches, mash_mode, max_dist));

            output_to!(distances, matches);
            return;
        }

        let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
        let mut filename_iter = filenames.iter();
        let filename = filename_iter.next().ok_or("At least one filename must be specified").unwrap();

        let first_sketch = open_mash_file(filename, matches, None).unwrap();

        // find the query sketches in the first file
        let mut query_sketches = Vec::new();
        if matches.is_present("queries") {
            let query_names: HashSet<String> = matches.values_of("queries").unwrap().map(|s| s.to_string()).collect();

            for sketch in first_sketch.sketches.iter() {
                if query_names.contains(&sketch.name) {
                    query_sketches.push(sketch.clone());
                }
            }
        } else {
            if first_sketch.sketches.len() < 1 {
                panic!("No sketches in first file!");
            }
            query_sketches.push(first_sketch.sketches[0].clone());
        }

        let mut distances = calc_sketch_distances(&query_sketches, &first_sketch.sketches, mash_mode, max_dist);
        for filename in filename_iter {
            let sketches = open_mash_file(filename, matches, Some(&first_sketch)).unwrap();
            distances.extend(calc_sketch_distances(&query_sketches, &sketches.sketches, mash_mode, max_dist));
        }
        output_to!(distances, matches);
    } else if let Some(matches) = matches.subcommand_matches("hist") {
        let mut hist_map: HashMap<String, Vec<u64>> = HashMap::new();
        parse_mash_files(matches, |multisketch, _| {
            for sketch in multisketch.sketches.iter() {
                hist_map.insert(sketch.name.to_string(), hist(&sketch.get_kmers().unwrap()));
            }
        }).unwrap();

        output_to!(hist_map, matches);
    } else if let Some(matches) = matches.subcommand_matches("info") {
        parse_mash_files(matches, |multisketch, _| {
            for sketch in multisketch.sketches.iter() {
                print!("{}", &sketch.name);
                if let Some(l) = sketch.seqLength {
                    println!(" (from {}bp)", l);
                } else {
                    println!("");
                }
                if let Some(kmers) = sketch.get_kmers() {
                    if let Ok(c) = cardinality(&kmers) {
                        println!("  Estimated # of Unique Kmers: {}", c);
                    }

                    let histogram = hist(&kmers);
                    let mean = histogram.iter().enumerate()
                        .map(|(i, v)| (i as f32 * *v as f32, *v as f32))
                        .fold((0f32, 0f32), |e, s| (e.0 + s.0, e.1 + s.1));
                    println!("  Estimated Average Depth: {}x", mean.0 / mean.1);

                    let mut total_gc: u64 = 0;
                    for kmer in &kmers {
                        total_gc += kmer.kmer.iter().map(|b| {
                            match *b {
                                b'G' | b'g' | b'C' | b'c' => kmer.count as u64,
                                _ => 0,
                            }
                        }).sum();
                    }
                    let total_bases = mean.0 * kmers[0].kmer.len() as f32;
                    println!("  Estimated % GC: {}%", 100f32 * total_gc as f32 / total_bases);
                }
            }
        }).unwrap();
    }
}


pub fn mash_files(filenames: Vec<&str>, n_hashes: usize, final_size: usize, kmer_length: u8, filters: &mut FilterParams, no_strict: bool, seed: u64) -> Result<JSONMultiSketch, String> {
    let mut sketches = Vec::with_capacity(filenames.len());
    for filename in &filenames {
        let mut seq_len = 0u64;
        let path = Path::new(filename);
        let mut minhash = match filters.filter_on {
            Some(true) | None => MinHashKmers::new(n_hashes, seed),
            Some(false) => MinHashKmers::new(final_size, seed),
        };
        fastx_cli(path.to_str().ok_or("Couldn't make path into string")?, |seq_type| {
            // disable filtering for FASTA files unless it was explicitly specified
            if let None = filters.filter_on {
                filters.filter_on = match seq_type {
                    "FASTA" => Some(false),
                    "FASTQ" => Some(true),
                    _ => panic!("Unknown sequence type"),
                };
            }
        }, |seq| {
            seq_len += seq.seq.len() as u64;
            for (_, kmer, is_rev_complement) in seq.normalize(false).kmers(kmer_length, true) {
                let rc_count = match is_rev_complement {
                    true => 1u8,
                    false => 0u8,
                };
                minhash.push(kmer, rc_count);
            }
        }).map_err(|e| e.to_string())?;

        let hashes = minhash.into_vec();
        let (mut filtered_hashes, filter_stats) = filter_sketch(&hashes, &filters);
        filtered_hashes.truncate(final_size);
        if !no_strict && filtered_hashes.len() < final_size {
            return Err(format!("{} had too few kmers ({}) to sketch", filename, filtered_hashes.len()));
        }

        // directory should be clipped from filename
        let basename = path.file_name().ok_or("Couldn't get filename from path")?;
        let sketch = JSONSketch::new(basename.to_str().ok_or("Couldn't make filename into string")?, seq_len, filtered_hashes, &filter_stats);
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

fn parse_all_mash_files(matches: &ArgMatches) -> Result<JSONMultiSketch, String> {
    let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
    let mut filename_iter = filenames.iter();
    let filename = filename_iter.next().ok_or("At least one filename must be specified")?;

    let mut sketches = open_mash_file(filename, matches, None)?;
    for filename in filename_iter {
        let sketch = open_mash_file(filename, matches, Some(&sketches))?;
        sketches.sketches.extend_from_slice(&sketch.sketches);
    }
    Ok(sketches)
}

fn parse_mash_files<P>(matches: &ArgMatches, ref mut parser: P) -> Result<(), String> where P: FnMut(&JSONMultiSketch, &str) -> () {
    let filenames: Vec<_> = matches.values_of("INPUT").unwrap().collect();
    let mut filename_iter = filenames.iter();
    let filename = filename_iter.next().ok_or("At least one filename must be specified")?;

    let first_sketch = open_mash_file(filename, matches, None)?;
    parser(&first_sketch, filename);
    for filename in filename_iter {
        let sketch = open_mash_file(filename, matches, Some(&first_sketch))?;
        parser(&sketch, filename);
    }
    Ok(())
}

fn calc_sketch_distances(ref_sketches: &[JSONSketch], query_sketches: &[JSONSketch], mash_mode: bool, max_distance: f64) -> Vec<SketchDistance> {
    let mut distances = Vec::new();
    for ref_sketch in ref_sketches.iter() {
        let rsketch = &ref_sketch.get_kmers().unwrap();
        for query_sketch in query_sketches.iter() {
            if query_sketch == ref_sketch {
                continue;
            }
            let qsketch = &query_sketch.get_kmers().unwrap();
            let distance = distance(&qsketch, &rsketch, &query_sketch.name, &ref_sketch.name, mash_mode).unwrap();
            if distance.mashDistance <= max_distance {
                distances.push(distance);
            }
        }
    }
    distances
}

fn open_mash_file(filename: &str, matches: &ArgMatches, default_sketch: Option<&JSONMultiSketch>) -> Result<JSONMultiSketch, String> {
    let final_sketch_size = matches.value_of("n_hashes").unwrap().parse::<usize>().map_err(|_| "n_hashes must be an integer")?;
    let seed = matches.value_of("seed").unwrap().parse::<u64>().map_err(|_| "seed must be an integer")?;
    let kmer_length = matches.value_of("kmer_length").unwrap().parse::<u8>().map_err(|_| "kmer_length must be an integer < 256")?;

    let no_strict = matches.is_present("no_strict");
    let filter_on = match (matches.is_present("filter"), matches.is_present("no_filter")) {
        (true, true) => panic!("Can't have both filtering and no filtering!"),
        (true, false) => Some(true),
        (false, true) => Some(false),
        (false, false) => None,
    };
    let min_abun_filter = match matches.occurrences_of("min_abun_filter") > 0 {
        true => Some(matches.value_of("min_abun_filter").unwrap().parse::<u16>().map_err(|_| {
            format!("min_abun_filter must be a number greater than or equal to 0")
        })?),
        false => None,
    };
    let max_abun_filter = match matches.occurrences_of("max_abun_filter") > 0 {
        true => Some(matches.value_of("max_abun_filter").unwrap().parse::<u16>().map_err(|_| {
            format!("max_abun_filter must be a number greater than or equal to 0")
        })?),
        false => None,
    };
    let err_filter = matches.value_of("err_filter").unwrap().parse::<f32>().map_err(|_| {
        format!("err-filter must be a number")
    }).and_then(|r| {
        if 0f32 <= r && r <= 100f32 / kmer_length as f32 {
            return Ok(kmer_length as f32 * r / 100f32);
        }
        Err(format!("err-filter must be a percent between 0 and {}", 100f32 / kmer_length as f32))
    })?;
    let strand_filter = matches.value_of("strand_filter").unwrap().parse::<f32>().map_err(|_| {
        format!("strand-filter must be a number")
    }).and_then(|r| {
        if 0f32 <= r && r <= 1f32 {
            return Ok(r);
        }
        Err(format!("strand-filter must be a ratio between 0 and 1"))
    })?;

    // note: is_present returns true while occurrences_of correctly is 0
    let mut filters = FilterParams {
        filter_on: filter_on,
        abun_filter: (min_abun_filter, max_abun_filter),
        err_filter: err_filter,
        strand_filter: strand_filter,
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
    let mut json: JSONMultiSketch = serde_json::from_reader(file).map_err(|_|
        format!("Error parsing {}", &filename)
    )?;

    // if filtering is explicitly set, re-filter the hashes
    if filters.filter_on == Some(true) {
        for sketch in &mut json.sketches {
            sketch.apply_filtering(&filters);
        }
    }

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
