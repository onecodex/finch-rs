use std::mem::size_of;
use std::str::FromStr;

use clap::{App, Arg, ArgMatches};
use failure::{bail, format_err};

use crate::filtering::FilterParams;
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
        get_float_arg::<f32>(matches, "err_filter", 100f64 / f64::from(kmer_length))?;
    err_filter *= f32::from(kmer_length) / 100f32;
    let strand_filter = get_float_arg::<f32>(matches, "strand_filter", 1f64)?;

    Ok(FilterParams {
        filter_on,
        abun_filter: (min_abun_filter, max_abun_filter),
        err_filter,
        strand_filter,
    })
}
