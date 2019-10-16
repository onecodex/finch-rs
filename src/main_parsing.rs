use std::mem::size_of;
use std::str::FromStr;

use clap::{App, Arg, ArgMatches};
use failure::{bail, format_err};

use crate::filtering::FilterParams;
use crate::serialization::MultiSketch;
use crate::sketch_schemes::SketchParams;
use crate::Result;

pub fn get_int_arg<T: FromStr>(matches: &ArgMatches, key: &str) -> Result<T> {
    let display_key = key.replace("_", "-");
    Ok(matches
        .value_of(key)
        .ok_or_else(|| format_err!("Bad {}", display_key))?
        .parse::<T>()
        .map_err(|_| {
            if size_of::<T>() == 1 {
                format_err!("{} must be a positive integer (<256)", display_key)
            } else {
                format_err!("{} must be a positive integer", display_key)
            }
        })?)
}

pub fn get_float_arg<T: Copy + FromStr + Into<f64>>(
    matches: &ArgMatches,
    key: &str,
    limit: f64,
) -> Result<T> {
    let display_key = key.replace("_", "-");
    Ok(matches
        .value_of(key)
        .ok_or_else(|| format_err!("Bad {}", display_key))?
        .parse::<T>()
        .map_err(|_| format_err!("{} must be a number", display_key))
        .and_then(|r| {
            if 0f64 <= r.into() && r.into() <= limit {
                return Ok(r);
            }
            bail!("{} must be between 0 and {}", display_key, limit)
        })?)
}

pub fn add_filter_options<'a, 'b>(app: App<'a, 'b>) -> App<'a, 'b> {
    app.arg(Arg::with_name("no_filter")
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
}

pub fn parse_filter_options(matches: &ArgMatches, kmer_length: u8) -> Result<FilterParams> {
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
        Some(get_int_arg::<u64>(matches, "min_abun_filter")?)
    } else {
        None
    };

    let max_abun_filter = if matches.occurrences_of("max_abun_filter") > 0 {
        Some(get_int_arg::<u64>(matches, "max_abun_filter")?)
    } else {
        None
    };

    let mut err_filter =
        get_float_arg::<f64>(matches, "err_filter", 100f64 / f64::from(kmer_length))?;
    err_filter *= f64::from(kmer_length) / 100f64;

    let strand_filter = get_float_arg::<f64>(matches, "strand_filter", 1f64)?;

    Ok(FilterParams {
        filter_on,
        abun_filter: (min_abun_filter, max_abun_filter),
        err_filter,
        strand_filter,
    })
}

pub fn add_sketch_options<'a, 'b>(app: App<'a, 'b>) -> App<'a, 'b> {
    // note we're defining groups for the arguments depending on which
    // sketch_type they're used with, but clap doesn't allow us to flag
    // argument conflicts/requirements based off another argument's value
    app.arg(Arg::with_name("sketch_type")
         .short("s")
         .long("sketch-type")
         .takes_value(true)
         .possible_values(&["mash", "scaled", "none"])
         .default_value("mash")
         .help("What type of sketching to perform"))
    .arg(Arg::with_name("kmer_length")
         .short("k")
         .long("kmer-length")
         .takes_value(true)
         .default_value_if("sketch_type", Some("none"), "4")
         .default_value("21")
         .help("Length of kmers to use"))
    .arg(Arg::with_name("n_hashes")
         .short("n")
         .long("n-hashes")
         .takes_value(true)
         // .groups(&["mash", "scaled"])
         .default_value("1000")
         .help("How many kmers/hashes to store [`sketch-type=mash` and `sketch-type=scaled`]"))
    .arg(Arg::with_name("scale")
         .long("scale")
         .takes_value(true)
         .default_value("0.001")
         // .group("scaled")
         .help("Sketch scaling factor [`sketch-type=scaled` only]"))
    .arg(Arg::with_name("seed")
         .long("seed")
         .takes_value(true)
         // .groups(&["mash", "scaled"])
         .default_value("0")
         .help("Seed murmurhash with this value [`sketch-type=mash` and `sketch-type=scaled`]"))
    .arg(Arg::with_name("oversketch")
         .long("oversketch")
         .takes_value(true)
         // .group("mash")
         .default_value("200")
         .help("The amount of extra sketching to do before filtering. This is only a safety to allow sketching e.g. high-coverage files with lots of error-generated uniquemers and should not change the final sketch [`sketch-type=mash` only]"))
    .arg(Arg::with_name("no_strict")
         .short("N")
         .long("no-strict")
         // .group("mash")
         .help("Allow sketching files with fewer kmers than `n_hashes` [`sketch-type=mash` only]"))
}

pub fn parse_sketch_options(
    matches: &ArgMatches,
    kmer_length: u8,
    filters_enabled: Option<bool>,
) -> Result<SketchParams> {
    Ok(match matches.value_of("sketch_type").unwrap_or("mash") {
        "mash" => {
            if matches.occurrences_of("scale") != 0 {
                bail!("`scale` can not be specified for `mash` sketch types")
            }
            let final_size: usize = get_int_arg(matches, "n_hashes")?;
            let oversketch: usize = get_int_arg(matches, "oversketch")?;
            let sketch_size = final_size * oversketch;

            let kmers_to_sketch = match filters_enabled {
                Some(true) | None => sketch_size,
                Some(false) => final_size,
            };

            SketchParams::Mash {
                kmers_to_sketch,
                final_size,
                no_strict: matches.is_present("no_strict"),
                kmer_length,
                hash_seed: get_int_arg(matches, "seed")?,
            }
        }
        "scaled" => {
            if matches.occurrences_of("oversketch") != 0 {
                bail!("`oversketch` can not be specified for `scaled` sketch types")
            }
            if matches.occurrences_of("no_strict") != 0 {
                bail!("`no_strict` can not be specified for `scaled` sketch types")
            }
            let kmers_to_sketch: usize = get_int_arg(matches, "n_hashes")?;
            let scale: f64 = get_float_arg(matches, "scale", 1.)?;
            SketchParams::Scaled {
                kmers_to_sketch,
                kmer_length,
                scale,
                hash_seed: get_int_arg(matches, "seed")?,
            }
        }
        "none" => {
            if matches.occurrences_of("n_hashes") != 0 {
                bail!("`n_hashes` can not be specified for `none` sketch types")
            }
            if matches.occurrences_of("seed") != 0 {
                bail!("`seed` can not be specified for `none` sketch types")
            }
            if matches.occurrences_of("oversketch") != 0 {
                bail!("`oversketch` can not be specified for `none` sketch types")
            }
            if matches.occurrences_of("no_strict") != 0 {
                bail!("`no_strict` can not be specified for `none` sketch types")
            }
            if matches.occurrences_of("scale") != 0 {
                bail!("`scale` can not be specified for `none` sketch types")
            }
            SketchParams::AllCounts { kmer_length }
        }
        _ => panic!("A unknown sketch type was selected"),
    })
}

pub fn update_sketch_params(
    matches: &ArgMatches,
    sketch_params: &mut SketchParams,
    multisketch: &MultiSketch,
    name: &str,
) -> Result<()> {
    // FIXME: check that the sketching type is concordant?

    // if arguments weren't provided use the ones from the multisketch
    match sketch_params {
        SketchParams::Mash {
            final_size,
            kmer_length,
            hash_seed,
            ..
        } => {
            if matches.occurrences_of("n_hashes") == 0 {
                *final_size = multisketch.sketch_size as usize;
            }
            if matches.occurrences_of("kmer_length") == 0 {
                *kmer_length = multisketch.kmer;
            } else if *kmer_length != multisketch.kmer {
                bail!(
                    "Specified kmer length {} does not match {} from sketch {}",
                    kmer_length,
                    multisketch.kmer,
                    name
                );
            }
            if matches.occurrences_of("seed") == 0 {
                *hash_seed = multisketch.hash_seed;
            } else if *hash_seed != multisketch.hash_seed {
                bail!(
                    "Specified hash seed {} does not match {} from sketch {}",
                    hash_seed,
                    multisketch.hash_seed,
                    name
                );
            }
        }
        _ => {}
    }
    Ok(())
}
