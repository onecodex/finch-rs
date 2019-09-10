#[macro_use]
extern crate clap;
#[macro_use]
extern crate failure;
extern crate finch;
extern crate serde_json;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{stderr, stdout, BufReader, Write};
use std::process::exit;

use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use memmap::MmapOptions;

use finch::distance::distance;
use finch::filtering::FilterParams;
use finch::serialization::{
    read_mash_file, write_mash_file, MultiSketch, Sketch, SketchDistance, FINCH_EXT, MASH_EXT,
};
use finch::statistics::{cardinality, hist};
use finch::{mash_files, Result};

macro_rules! add_output_options {
    ($cmd:ident) => {
        $cmd = $cmd
            .arg(
                Arg::with_name("output_file")
                    .short("o")
                    .long("output")
                    .help("Output to this file")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("std_out")
                    .short("O")
                    .long("std-out")
                    .help("Output to stdout ('print to terminal')")
                    .conflicts_with("output_file"),
            );
    };
}

fn output_to<F>(output_fn: F, matches: &ArgMatches, extension: &str) -> Result<()>
where
    F: Fn(&mut dyn Write) -> Result<()>,
{
    let output = matches.value_of("output_file");
    match output {
        None => {
            let mut out = stdout();
            output_fn(&mut out)?;
        }
        Some(o) => {
            // if the filename doesn't have the right extension
            // add it on
            let filename = String::from(o);
            let out_filename = if filename.ends_with(extension) {
                filename
            } else {
                filename + extension
            };

            let mut out = File::create(&out_filename)
                .map_err(|_| format_err!("Could not create {}", out_filename))?;
            output_fn(&mut out)?;
        }
    };
    Ok(())
}

macro_rules! add_kmer_options {
    ($cmd:ident) => {
        $cmd = $cmd.arg(Arg::with_name("n_hashes")
             .short("n")
             .long("n-hashes")
             .help("How many kmers/hashes to store")
             .takes_value(true)
             .default_value("1000"))
        .arg(Arg::with_name("scaled")
             .short("s")
             .long("scaled")
             .help("Approximate percentage of total kmers/hashes to store")
             .takes_value(true))
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
    // see https://github.com/rust-lang-nursery/failure/issues/76
    if let Err(err) = run() {
        let mut serr = stderr();
        let mut causes = err.iter_chain();
        writeln!(
            serr,
            "Error: {}",
            causes.next().expect("`causes` to at least contain `err`")
        )
        .expect("unable to write error to stderr");
        for cause in causes {
            writeln!(serr, "Caused by: {}", cause).expect("unable to write error to stderr");
        }
        // The following assumes an `Error`, use `if let Some(backtrace) ...` for a `Fail`
        writeln!(serr, "{:?}", err.backtrace()).expect("unable to write error to stderr");
        exit(1);
    }
}

fn run() -> Result<()> {
    let mut sketch_command = SubCommand::with_name("sketch")
        .about("Create sketches from FASTA/Q file(s)")
        .arg(
            Arg::with_name("INPUT")
                .help("The file(s) to sketch")
                .multiple(true)
                .required(true),
        );
    if cfg!(feature = "mash_format") {
        sketch_command = sketch_command.arg(
            Arg::with_name("binary_format")
                .short("b")
                .long("binary-format")
                .help("Outputs sketch in a binary format compatible with `mash`"),
        );
    }
    add_output_options!(sketch_command);
    add_kmer_options!(sketch_command);

    let mut dist_command = SubCommand::with_name("dist")
        .about("Compute distances between sketches")
        .arg(
            Arg::with_name("INPUT")
                .help("Sketchfile(s) to make comparisons for")
                .multiple(true)
                .required(true),
        )
        .arg(
            Arg::with_name("pairwise")
                .short("p")
                .long("pairwise")
                .conflicts_with("queries")
                .help("Calculate distances between all sketches"),
        )
        .arg(
            Arg::with_name("queries")
                .short("q")
                .long("queries")
                .help("All distances are from these sketches (sketches must be in the first file)")
                .multiple(true)
                .conflicts_with("pairwise")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max_distance")
                .short("d")
                .long("max-dist")
                .help("Only report distances under this threshold")
                .default_value("1.0")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mash_mode")
                .short("m")
                .long("mash")
                .help("Calculate distances using the same algorithms as Mash"),
        );
    add_output_options!(dist_command);
    add_kmer_options!(dist_command);

    let mut hist_command = SubCommand::with_name("hist")
        .about("Display histograms of kmer abundances")
        .arg(
            Arg::with_name("INPUT")
                .help("Generate histograms from these file(s)")
                .multiple(true)
                .required(true),
        );
    add_output_options!(hist_command);
    add_kmer_options!(hist_command);

    let mut info_command = SubCommand::with_name("info")
        .about("Display basic statistics")
        .arg(
            Arg::with_name("INPUT")
                .help("Return stats on these file(s)")
                .multiple(true)
                .required(true),
        );
    add_output_options!(info_command);
    add_kmer_options!(info_command);

    let matches = App::new("finch")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Tool for working with genomic MinHash sketches")
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(sketch_command)
        .subcommand(dist_command)
        .subcommand(hist_command)
        .subcommand(info_command)
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("sketch") {
        let file_ext = if matches.is_present("binary_format") {
            MASH_EXT
        } else {
            FINCH_EXT
        };
        if matches.is_present("output_file") || matches.is_present("std_out") {
            let sketches = parse_all_mash_files(matches)?;
            output_to(
                |writer| {
                    if matches.is_present("binary_format") {
                        write_mash_file(writer, &sketches)?;
                    } else {
                        serde_json::to_writer(writer, &sketches)?;
                    }
                    Ok(())
                },
                &matches,
                &file_ext,
            )?;
        } else {
            // "sketch in place"
            parse_mash_files(matches, |multisketch, filename| {
                let out_filename = filename.to_string() + file_ext;
                let mut out = File::create(&out_filename)
                    .map_err(|_| format_err!("Could not open {}", out_filename))
                    .unwrap();
                if matches.is_present("binary_format") {
                    write_mash_file(&mut out, &multisketch).unwrap();
                } else {
                    serde_json::to_writer(&mut out, &multisketch).unwrap();
                }
            })?;
        }
    } else if let Some(matches) = matches.subcommand_matches("dist") {
        let mash_mode = matches.is_present("mash_mode");

        let max_dist = matches
            .value_of("max_distance")
            .ok_or_else(|| format_err!("Bad max-distance"))?
            .parse::<f64>()
            .map_err(|_| format_err!("max-distance must be a number"))
            .and_then(|r| {
                if 0f64 <= r && r <= 1f64 {
                    return Ok(r);
                }
                bail!("max-distance must be 0 and 1")
            })?;

        let all_sketches = parse_all_mash_files(matches)?;

        let mut query_sketches = Vec::new();
        if matches.is_present("pairwise") {
            for sketch in &all_sketches.sketches {
                query_sketches.push(sketch);
            }
        } else if matches.is_present("queries") {
            let query_names: HashSet<String> = matches
                .values_of("queries")
                .ok_or_else(|| format_err!("Bad queries"))?
                .map(|s| s.to_string())
                .collect();

            for sketch in all_sketches.sketches.iter() {
                if query_names.contains(&sketch.name) {
                    query_sketches.push(sketch);
                }
            }
        } else {
            if all_sketches.sketches.is_empty() {
                bail!("No sketches present!");
            }
            query_sketches.push(all_sketches.sketches.last().unwrap());
        }

        let distances =
            calc_sketch_distances(&query_sketches, &all_sketches.sketches, mash_mode, max_dist);

        output_to(
            |writer| {
                serde_json::to_writer(writer, &distances)
                    .map_err(|_| format_err!("Could not serialize JSON to file"))?;
                Ok(())
            },
            &matches,
            ".json",
        )?;
    } else if let Some(matches) = matches.subcommand_matches("hist") {
        let mut hist_map: HashMap<String, Vec<u64>> = HashMap::new();
        parse_mash_files(matches, |multisketch, _| {
            for sketch in multisketch.sketches.iter() {
                hist_map.insert(sketch.name.to_string(), hist(&sketch.hashes));
            }
        })?;

        output_to(
            |writer| {
                serde_json::to_writer(writer, &hist_map)
                    .map_err(|_| format_err!("Could not serialize JSON to file"))?;
                Ok(())
            },
            &matches,
            ".json",
        )?;
    } else if let Some(matches) = matches.subcommand_matches("info") {
        // TODO: this should probably output JSON
        parse_mash_files(matches, |multisketch, _| {
            for sketch in multisketch.sketches.iter() {
                print!("{}", &sketch.name);
                if let Some(l) = sketch.seqLength {
                    println!(" (from {}bp)", l);
                } else {
                    println!();
                }
                let kmers = &sketch.hashes;
                if let Ok(c) = cardinality(kmers) {
                    println!("  Estimated # of Unique Kmers: {}", c);
                }

                let histogram = hist(kmers);
                let mean = histogram
                    .iter()
                    .enumerate()
                    .map(|(i, v)| ((i as f32 + 1f32) * *v as f32, *v as f32))
                    .fold((0f32, 0f32), |e, s| (e.0 + s.0, e.1 + s.1));
                println!("  Estimated Average Depth: {}x", mean.0 / mean.1);

                let mut total_gc: u64 = 0;
                for kmer in kmers {
                    total_gc += kmer
                        .kmer
                        .iter()
                        .map(|b| match *b {
                            b'G' | b'g' | b'C' | b'c' => u64::from(kmer.count),
                            _ => 0,
                        })
                        .sum::<u64>();
                }
                let total_bases = if kmers.is_empty() {
                    0f32
                } else {
                    mean.0 * kmers[0].kmer.len() as f32
                };
                println!(
                    "  Estimated % GC: {}%",
                    100f32 * total_gc as f32 / total_bases
                );
            }
        })?;
    }
    Ok(())
}

fn parse_all_mash_files(matches: &ArgMatches) -> Result<MultiSketch> {
    let filenames: Vec<_> = matches
        .values_of("INPUT")
        .ok_or_else(|| format_err!("Bad INPUT"))?
        .collect();
    let mut filename_iter = filenames.iter();
    let filename = filename_iter
        .next()
        .ok_or_else(|| format_err!("At least one filename must be specified"))?;

    let mut sketches = open_mash_file(filename, matches, None)?;
    for filename in filename_iter {
        let sketch = open_mash_file(filename, matches, Some(&sketches))?;
        sketches.sketches.extend_from_slice(&sketch.sketches);
    }
    Ok(sketches)
}

fn parse_mash_files<P>(matches: &ArgMatches, mut parser: P) -> Result<()>
where
    P: FnMut(&MultiSketch, &str) -> (),
{
    let filenames: Vec<_> = matches
        .values_of("INPUT")
        .ok_or_else(|| format_err!("Bad INPUT"))?
        .collect();
    let mut filename_iter = filenames.iter();
    let filename = filename_iter
        .next()
        .ok_or_else(|| format_err!("At least one filename must be specified"))?;

    let first_sketch = open_mash_file(filename, matches, None)?;
    parser(&first_sketch, filename);
    for filename in filename_iter {
        let sketch = open_mash_file(filename, matches, Some(&first_sketch))?;
        parser(&sketch, filename);
    }
    Ok(())
}

fn calc_sketch_distances(
    query_sketches: &[&Sketch],
    ref_sketches: &[Sketch],
    mash_mode: bool,
    max_distance: f64,
) -> Vec<SketchDistance> {
    let mut distances = Vec::new();
    for ref_sketch in ref_sketches.iter() {
        let rsketch = &ref_sketch.hashes;
        for query_sketch in query_sketches.iter() {
            if query_sketch == &ref_sketch {
                continue;
            }
            let qsketch = &query_sketch.hashes;
            let distance = distance(
                &qsketch,
                &rsketch,
                &query_sketch.name,
                &ref_sketch.name,
                mash_mode,
            )
            .unwrap();
            if distance.mashDistance <= max_distance {
                distances.push(distance);
            }
        }
    }
    distances
}

fn open_mash_file(
    filename: &str,
    matches: &ArgMatches,
    default_sketch: Option<&MultiSketch>,
) -> Result<MultiSketch> {
    let final_sketch_size = matches
        .value_of("n_hashes")
        .ok_or_else(|| format_err!("Bad n_hashes"))?
        .parse::<usize>()
        .map_err(|_| format_err!("n_hashes must be an integer"))?;
    let seed = matches
        .value_of("seed")
        .ok_or_else(|| format_err!("Bad seed"))?
        .parse::<u64>()
        .map_err(|_| format_err!("seed must be an integer"))?;
    let kmer_length = matches
        .value_of("kmer_length")
        .ok_or_else(|| format_err!("Bad kmer_length"))?
        .parse::<u8>()
        .map_err(|_| format_err!("kmer_length must be an integer < 256"))?;

    let scaled = if matches.occurrences_of("scaled") > 0 {
        Some(
            matches
                .value_of("scaled")
                .ok_or_else(|| format_err!("Bad scaled value"))?
                .parse::<f64>()
                .map_err(|_| format_err!("scaled must be in range (0, 1]"))?,
        )
    } else {
        None
    };

    let no_strict = matches.is_present("no_strict");
    let filter_on = match (
        matches.is_present("filter"),
        matches.is_present("no_filter"),
    ) {
        (true, true) => panic!("Can't have both filtering and no filtering!"),
        (true, false) => Some(true),
        (false, true) => Some(false),
        (false, false) => None,
    };
    let min_abun_filter = if matches.occurrences_of("min_abun_filter") > 0 {
        Some(
            matches
                .value_of("min_abun_filter")
                .ok_or_else(|| format_err!("Bad min_abun_filter"))?
                .parse::<u16>()
                .map_err(|_| {
                    format_err!("min_abun_filter must be a number greater than or equal to 0")
                })?,
        )
    } else {
        None
    };
    let max_abun_filter = if matches.occurrences_of("max_abun_filter") > 0 {
        Some(
            matches
                .value_of("max_abun_filter")
                .ok_or_else(|| format_err!("Bad max_abun_filter"))?
                .parse::<u16>()
                .map_err(|_| {
                    format_err!("max_abun_filter must be a number greater than or equal to 0")
                })?,
        )
    } else {
        None
    };
    let err_filter = matches
        .value_of("err_filter")
        .ok_or_else(|| format_err!("Bad err_filter"))?
        .parse::<f32>()
        .map_err(|_| format_err!("err-filter must be a number"))
        .and_then(|r| {
            if 0f32 <= r && r <= 100f32 / f32::from(kmer_length) {
                return Ok(f32::from(kmer_length) * r / 100f32);
            }
            bail!(
                "err-filter must be a percent between 0 and {}",
                100f32 / f32::from(kmer_length)
            )
        })?;
    let strand_filter = matches
        .value_of("strand_filter")
        .ok_or_else(|| format_err!("Bad strand_filter"))?
        .parse::<f32>()
        .map_err(|_| format_err!("strand-filter must be a number"))
        .and_then(|r| {
            if 0f32 <= r && r <= 1f32 {
                return Ok(r);
            }
            bail!("strand-filter must be a ratio between 0 and 1")
        })?;

    // note: is_present returns true while occurrences_of correctly is 0
    let mut filters = FilterParams {
        filter_on,
        abun_filter: (min_abun_filter, max_abun_filter),
        err_filter,
        strand_filter,
    };

    let oversketch = matches
        .value_of("oversketch")
        .ok_or_else(|| format_err!("Bad oversketch"))?
        .parse::<usize>()
        .map_err(|_| format_err!("bad value for oversketch"))?;
    let sketch_size = final_sketch_size * oversketch;

    // if the file isn't a sketch file, we try to sketch it and pass the sketches back
    if !filename.ends_with(".json")
        && !filename.ends_with(FINCH_EXT)
        && !filename.ends_with(MASH_EXT)
    {
        return match default_sketch {
            Some(s) => mash_files(
                &[filename],
                sketch_size / final_sketch_size * s.sketchSize as usize,
                s.sketchSize as usize,
                s.kmer,
                &mut filters,
                no_strict,
                s.hashSeed,
                scaled,
            ),
            None => mash_files(
                &[filename],
                sketch_size,
                final_sketch_size,
                kmer_length,
                &mut filters,
                no_strict,
                seed,
                scaled,
            ),
        };
    }

    // otherwise we just open the file and return the sketches
    let file = File::open(filename).map_err(|_| format_err!("Error opening {}", &filename))?;
    let mut json: MultiSketch = if filename.ends_with(MASH_EXT) {
        let mut buf_reader = BufReader::new(file);
        read_mash_file(&mut buf_reader)?
    } else {
        let mapped = unsafe { MmapOptions::new().map(&file)? };
        serde_json::from_slice(&mapped).map_err(|_| format_err!("Error parsing {}", &filename))?
    };

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
            bail!(
                "{} has a different kmer length ({}) from others ({})",
                filename,
                s.kmer,
                json.kmer
            );
        } else if s.hashType != json.hashType {
            bail!(
                "{} used a different hash ({}) from others ({})",
                filename,
                s.hashType,
                json.hashType
            );
        } else if s.hashSeed != json.hashSeed {
            bail!(
                "{} had a different hash seed ({}) from others ({})",
                filename,
                s.hashSeed,
                json.hashSeed
            );
        } else if s.hashBits != json.hashBits {
            bail!(
                "{} used a different length hash ({}) from others ({})",
                filename,
                s.hashBits,
                json.hashBits
            );
        }
    }
    Ok(json)
}
