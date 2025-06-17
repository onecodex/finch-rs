use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{stdout, Write};

use anyhow::{anyhow, bail, Context, Result};
use clap::ArgMatches;

use crate::cli::{get_float_arg, get_int_arg, parse_filter_options, parse_sketch_options};
use finch::distance::distance;
use finch::serialization::{
    write_finch_file, write_mash_file, MultiSketch, Sketch, SketchDistance, FINCH_BIN_EXT,
    FINCH_EXT, MASH_EXT,
};
use finch::sketch_schemes::SketchParams;
use finch::statistics::{cardinality, hist};
use finch::{format_err, open_sketch_file, sketch_files};
use std::mem::discriminant;

mod cli;

fn output_to<F>(output_fn: F, output: Option<&str>, extension: &str) -> Result<()>
where
    F: Fn(&mut dyn Write) -> Result<()>,
{
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
                .context(format!("unable to create '{}'", out_filename))?;
            output_fn(&mut out)?;
        }
    };
    Ok(())
}

fn run() -> Result<()> {
    let matches = cli::build_cli().get_matches();

    match matches.subcommand() {
        ("sketch", Some(matches)) => {
            let file_ext = if matches.is_present("binary_format") {
                FINCH_BIN_EXT
            } else if matches.is_present("mash_binary_format") {
                MASH_EXT
            } else {
                FINCH_EXT
            };

            if matches.is_present("output_file") || matches.is_present("std_out") {
                let sketches = parse_mash_files(matches)?;
                let output = matches.value_of("output_file");

                output_to(
                    |writer| {
                        if file_ext == FINCH_BIN_EXT {
                            write_finch_file(writer, &sketches)?;
                        } else if file_ext == MASH_EXT {
                            write_mash_file(writer, &sketches)?;
                        } else {
                            let multisketch = MultiSketch::from_sketches(&sketches)?;
                            serde_json::to_writer(writer, &multisketch)?;
                        }
                        Ok(())
                    },
                    output,
                    file_ext,
                )?;
            } else {
                // special case for "sketching in place"
                generate_sketch_files(matches, file_ext)?;
            }
        }
        ("dist", Some(matches)) => {
            let old_mode = matches.is_present("old_dist_mode");

            let max_dist = get_float_arg(matches, "max_distance", 1f64)?;
            let all_sketches = parse_mash_files(matches)?;

            let mut query_sketches = Vec::new();
            if matches.is_present("pairwise") {
                for sketch in &all_sketches {
                    query_sketches.push(sketch);
                }
            } else if matches.is_present("queries") {
                let query_names: HashSet<String> = matches
                    .values_of("queries")
                    .unwrap() // we already know it's present
                    .map(|s| s.to_string())
                    .collect();

                for sketch in &all_sketches {
                    if query_names.contains(&sketch.name) {
                        query_sketches.push(sketch);
                    }
                }
            } else {
                if all_sketches.is_empty() {
                    bail!("No sketches present!");
                }
                query_sketches.push(all_sketches.first().unwrap());
            }

            let distances =
                calc_sketch_distances(&query_sketches, &all_sketches, old_mode, max_dist);

            output_to(
                |writer| {
                    serde_json::to_writer(writer, &distances)
                        .map_err(|_| anyhow!("Could not serialize JSON to file"))?;
                    Ok(())
                },
                matches.value_of("output_file"),
                ".json",
            )?;
        }
        ("hist", Some(matches)) => {
            let mut hist_map: HashMap<String, Vec<u64>> = HashMap::new();
            let multisketch = parse_mash_files(matches)?;

            for sketch in multisketch {
                hist_map.insert(sketch.name.to_string(), hist(&sketch.hashes));
            }

            output_to(
                |writer| {
                    serde_json::to_writer(writer, &hist_map)
                        .map_err(|_| anyhow!("Could not serialize JSON to file"))?;
                    Ok(())
                },
                matches.value_of("output_file"),
                ".json",
            )?;
        }
        ("info", Some(matches)) => {
            // TODO: this should probably output JSON
            let multisketch = parse_mash_files(matches)?;

            for sketch in multisketch {
                print!("{}", &sketch.name);
                println!(" (from {}bp)", sketch.seq_length);
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
        }
        other => bail!("Unknown subcommand: {:?}", other.0),
    };

    Ok(())
}

fn main() {
    if let Err(err) = run() {
        eprintln!("Error: {:?}", err);
        std::process::exit(1);
    }
}

fn generate_sketch_files(matches: &ArgMatches, file_ext: &str) -> Result<()> {
    let filenames: Vec<_> = matches
        .values_of("INPUT")
        .ok_or_else(|| format_err!("Bad INPUT"))?
        .collect();

    let kmer_length: u8 = get_int_arg(matches, "kmer_length")?;
    let filters = parse_filter_options(matches, kmer_length)?;
    let sketch_params = parse_sketch_options(matches, kmer_length, filters.filter_on)?;

    for filename in filenames {
        if filename.ends_with(".json")
            || filename.ends_with(FINCH_EXT)
            || filename.ends_with(FINCH_BIN_EXT)
            || filename.ends_with(MASH_EXT)
        {
            bail!("Filename {} is not a sequence file?", filename);
        }

        let sketches = sketch_files(&[filename], &sketch_params, &filters)?;

        let out_filename = filename.to_string() + file_ext;
        let mut out = File::create(&out_filename)
            .map_err(|_| format_err!("Could not open {}", out_filename))?;
        if matches.is_present("binary_format") {
            write_finch_file(&mut out, &sketches)?;
        } else if matches.is_present("mash_binary_format") {
            write_mash_file(&mut out, &sketches)?;
        } else {
            let multisketch: MultiSketch = MultiSketch::from_sketches(&sketches)?;
            serde_json::to_writer(&mut out, &multisketch)?;
        }
    }
    Ok(())
}

fn parse_mash_files(matches: &ArgMatches) -> Result<Vec<Sketch>> {
    let filenames: Vec<_> = matches
        .values_of("INPUT")
        .ok_or_else(|| anyhow!("Bad INPUT"))?
        .collect();

    let mut sketch_filenames = Vec::new();
    let mut seq_filenames = Vec::new();
    for filename in filenames {
        if filename.ends_with(".json")
            || filename.ends_with(FINCH_EXT)
            || filename.ends_with(FINCH_BIN_EXT)
            || filename.ends_with(MASH_EXT)
        {
            sketch_filenames.push(filename);
        } else {
            seq_filenames.push(filename);
        }
    }

    let kmer_length: u8 = get_int_arg(matches, "kmer_length")?;
    let mut filters = parse_filter_options(matches, kmer_length)?;
    let mut sketch_params = parse_sketch_options(matches, kmer_length, filters.filter_on)?;

    let mut filename_iter = sketch_filenames.iter();
    if let Some(first_filename) = filename_iter.next() {
        let mut sketches = open_sketch_file(first_filename)?;

        update_sketch_params(matches, &mut sketch_params, &sketches[0], first_filename)?;
        // we also have to handle updating filter options separately because
        // kmer_length changes how we calculate the `err_filter`
        if matches.occurrences_of("kmer_length") == 0 {
            filters = parse_filter_options(matches, sketch_params.k())?;
        }

        // now do the filtering for the first sketch file
        if filters.filter_on == Some(true) {
            for sketch in &mut sketches {
                filters.filter_sketch(sketch);
            }
        }

        // and then handle the rest of the sketch files
        for filename in filename_iter {
            let extra_sketches = open_sketch_file(filename)?;
            // check new sketches are compatible with original file
            for sketch in &extra_sketches {
                if let Some((name, v1, v2)) =
                    sketch_params.check_compatibility(&sketch.sketch_params)
                {
                    bail!(
                        "Sketch {} has {} {}, but working value is {}",
                        sketch.name,
                        name,
                        v2,
                        v1,
                    );
                }
            }
            sketches.extend(extra_sketches);
            if filters.filter_on == Some(true) {
                for sketch in &mut sketches {
                    filters.filter_sketch(sketch);
                }
            }
        }

        // now handle the sequences
        let extra_sketches = sketch_files(&seq_filenames, &sketch_params, &filters)?;
        sketches.extend(extra_sketches);
        Ok(sketches)
    } else {
        // now handle the sequences
        let sketches = sketch_files(&seq_filenames, &sketch_params, &filters)?;
        Ok(sketches)
    }
}

fn calc_sketch_distances(
    query_sketches: &[&Sketch],
    ref_sketches: &[Sketch],
    old_mode: bool,
    max_distance: f64,
) -> Vec<SketchDistance> {
    let mut distances = Vec::new();
    for ref_sketch in ref_sketches {
        for query_sketch in query_sketches {
            if query_sketch == &ref_sketch {
                continue;
            }
            let distance = distance(query_sketch, ref_sketch, old_mode).unwrap();
            if distance.mash_distance <= max_distance {
                distances.push(distance);
            }
        }
    }
    distances
}

pub fn update_sketch_params(
    matches: &ArgMatches,
    sketch_params: &mut SketchParams,
    sketch: &Sketch,
    name: &str,
) -> Result<()> {
    let new_sketch_params = &sketch.sketch_params;

    // check that the sketching type is the same; we may want to remove
    // this check at some point?
    if discriminant(sketch_params) != discriminant(new_sketch_params) {
        bail!("Sketch types are not the same")
    }

    // if arguments weren't provided use the ones from the multisketch
    match sketch_params {
        SketchParams::Mash {
            final_size,
            kmer_length,
            hash_seed,
            ..
        } => {
            if matches.occurrences_of("n_hashes") == 0 {
                *final_size = new_sketch_params.expected_size();
            }
            if matches.occurrences_of("kmer_length") == 0 {
                *kmer_length = new_sketch_params.k();
            } else if *kmer_length != new_sketch_params.k() {
                bail!(
                    "Specified kmer length {} does not match {} from sketch {}",
                    kmer_length,
                    new_sketch_params.k(),
                    name
                );
            }
            let (_, _, new_hash_seed, _) = new_sketch_params.hash_info();
            if matches.occurrences_of("seed") == 0 {
                *hash_seed = new_hash_seed;
            } else if *hash_seed != new_hash_seed {
                bail!(
                    "Specified hash seed {} does not match {} from sketch {}",
                    hash_seed,
                    new_hash_seed,
                    name
                );
            }
            // TODO: do we need to update any of the other sketch params?
        }
        SketchParams::Scaled {
            kmer_length,
            hash_seed,
            scale,
            ..
        } => {
            // TODO: these two are identical to above so some DRY might be good?
            if matches.occurrences_of("kmer_length") == 0 {
                *kmer_length = new_sketch_params.k();
            } else if *kmer_length != new_sketch_params.k() {
                bail!(
                    "Specified kmer length {} does not match {} from sketch {}",
                    kmer_length,
                    new_sketch_params.k(),
                    name
                );
            }
            let (_, _, new_hash_seed, new_scale) = new_sketch_params.hash_info();
            if matches.occurrences_of("seed") == 0 {
                *hash_seed = new_hash_seed;
            } else if *hash_seed != new_hash_seed {
                bail!(
                    "Specified hash seed {} does not match {} from sketch {}",
                    hash_seed,
                    new_hash_seed,
                    name
                );
            }

            if let Some(new_scale_num) = new_scale {
                if matches.occurrences_of("scale") == 0 {
                    *scale = new_scale_num;
                } else if (*scale - new_scale_num).abs() < f64::EPSILON {
                    // TODO: maybe this should have a slightly larger delta?
                    bail!(
                        "Specified scale {} does not match {} from sketch {}",
                        scale,
                        new_scale_num,
                        name
                    );
                }
            }
        }
        SketchParams::AllCounts { kmer_length } => {
            if matches.occurrences_of("kmer_length") == 0 {
                *kmer_length = new_sketch_params.k();
            } else if *kmer_length != new_sketch_params.k() {
                bail!(
                    "Specified kmer length {} does not match {} from sketch {}",
                    kmer_length,
                    new_sketch_params.k(),
                    name
                );
            }
        }
    }
    Ok(())
}
