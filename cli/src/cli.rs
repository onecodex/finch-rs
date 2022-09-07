use anyhow::{anyhow, bail, Result};
use clap::{crate_authors, crate_version, App, AppSettings, Arg, ArgMatches, SubCommand};
use finch::filtering::FilterParams;
use finch::sketch_schemes::SketchParams;
use std::str::FromStr;

pub fn build_cli() -> App<'static, 'static> {
    App::new("finch")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Tool for working with genomic MinHash sketches")
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(info_command())
        .subcommand(sketch_command())
        .subcommand(dist_command())
        .subcommand(hist_command())
}

fn info_command() -> App<'static, 'static> {
    let mut info_command = SubCommand::with_name("info")
        .about("Display basic statistics")
        .arg(
            Arg::with_name("INPUT")
                .help("Return stats on these file(s)")
                .multiple(true)
                .required(true),
        );
    info_command = add_filter_options(info_command);
    info_command = add_sketch_options(info_command);
    info_command
}

fn sketch_command() -> App<'static, 'static> {
    let mut sketch_command = SubCommand::with_name("sketch")
        .about("Create sketches from FASTA/Q file(s)")
        .arg(
            Arg::with_name("INPUT")
                .help("The file(s) to sketch")
                .multiple(true)
                .required(true),
        )
        .arg(
            Arg::with_name("binary_format")
                .short("b")
                .long("finch-binary-format")
                .help("Outputs sketch to a finch-native binary format"),
        )
        .arg(
            Arg::with_name("mash_binary_format")
                .short("B")
                .long("mash-binary-format")
                .conflicts_with("binary_format")
                .help("Outputs sketch in a binary format compatible with `mash`"),
        );
    sketch_command = add_output_options(sketch_command);
    sketch_command = add_filter_options(sketch_command);
    sketch_command = add_sketch_options(sketch_command);
    sketch_command
}

fn dist_command() -> App<'static, 'static> {
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
            Arg::with_name("old_dist_mode")
                .long("old-dist")
                .help("Calculate distances using the old containment-biased Finch mode"),
        );
    dist_command = add_output_options(dist_command);
    dist_command = add_filter_options(dist_command);
    dist_command = add_sketch_options(dist_command);
    dist_command
}

fn hist_command() -> App<'static, 'static> {
    let mut hist_command = SubCommand::with_name("hist")
        .about("Display histograms of kmer abundances")
        .arg(
            Arg::with_name("INPUT")
                .help("Generate histograms from these file(s)")
                .multiple(true)
                .required(true),
        );
    hist_command = add_output_options(hist_command);
    hist_command = add_filter_options(hist_command);
    hist_command = add_sketch_options(hist_command);
    hist_command
}

fn add_filter_options<'a, 'b>(app: App<'a, 'b>) -> App<'a, 'b> {
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
         .help("The assumed error rate (as a percentage) used to dynamically determine the minimum coverage threshold from the kmer count histogram. This threshold is then used in place of the min-abun-filter if it's more stringent.")
         .takes_value(true)
         .default_value("1"))
}

fn add_sketch_options<'a, 'b>(app: App<'a, 'b>) -> App<'a, 'b> {
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

fn add_output_options<'a, 'b>(app: App<'a, 'b>) -> App<'a, 'b> {
    app.arg(
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
    )
}

pub fn get_int_arg<T: FromStr>(matches: &ArgMatches, key: &str) -> Result<T> {
    let display_key = key.replace('_', "-");
    matches
        .value_of(key)
        .ok_or_else(|| anyhow!("Bad {}", display_key))?
        .parse::<T>()
        .map_err(|_| anyhow!("{} must be a positive integer", display_key))
}

pub fn get_float_arg(matches: &ArgMatches, key: &str, limit: f64) -> Result<f64> {
    let display_key = key.replace('_', "-");
    matches
        .value_of(key)
        .ok_or_else(|| anyhow!("Bad {}", display_key))?
        .parse::<f64>()
        .map_err(|_| anyhow!("{} must be a number", display_key))
        .and_then(|r| {
            if 0f64 <= r && r <= limit {
                return Ok(r);
            }
            bail!("{} must be between 0 and {}", display_key, limit)
        })
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
        Some(get_int_arg::<u32>(matches, "min_abun_filter")?)
    } else {
        None
    };

    let max_abun_filter = if matches.occurrences_of("max_abun_filter") > 0 {
        Some(get_int_arg::<u32>(matches, "max_abun_filter")?)
    } else {
        None
    };

    let mut err_filter = get_float_arg(matches, "err_filter", 100f64 / f64::from(kmer_length))?;
    err_filter *= f64::from(kmer_length) / 100f64;

    let strand_filter = get_float_arg(matches, "strand_filter", 1f64)?;

    Ok(FilterParams {
        filter_on,
        abun_filter: (min_abun_filter, max_abun_filter),
        err_filter,
        strand_filter,
    })
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
        _ => bail!("A unknown sketch type was selected"),
    })
}
