#![cfg_attr(feature = "python", feature(specialization))]

#[cfg(feature = "mash_format")]
extern crate capnp;
#[macro_use]
extern crate serde_derive;

use std::fs::File;
use std::io::{stdin, Read};
use std::path::Path;
use std::result::Result as StdResult;

use failure::{format_err, Error};
use needletail::formats::parse_sequence_reader;
use rayon::prelude::*;

use crate::filtering::FilterParams;
use crate::serialization::{MultiSketch, Sketch};
use crate::sketch_schemes::SketchParams;

pub mod distance;
pub mod filtering;
pub mod sketch_schemes;
// it would be nice if there was a `pub(in main)` or something for
// main_parsing so we don't import it for `lib` itself
pub mod main_parsing;
#[cfg(feature = "mash_format")]
mod mash_capnp;
#[cfg(feature = "python")]
pub mod python;
pub mod serialization;
pub mod statistics;

pub type Result<T> = StdResult<T, Error>;

pub fn sketch_files(
    filenames: &[&str],
    sketch_params: &SketchParams,
    filters: &FilterParams,
) -> Result<MultiSketch> {
    let sketches: Result<Vec<Sketch>> = filenames
        .par_iter()
        .map(|filename| {
            // open the file with a special case to handle stdin
            let sin = stdin();
            let reader: Box<dyn Read> = if filename == &"-" {
                Box::new(sin.lock())
            } else {
                Box::new(File::open(&Path::new(filename))?)
            };
            // sketch!
            Ok(sketch_stream(reader, filename, sketch_params, &filters)?)
        })
        .collect();
    let (hash_type, hash_bits, hash_seed) = sketch_params.hash_info();
    Ok(MultiSketch {
        kmer: sketch_params.k(),
        alphabet: String::from("ACGT"),
        preserve_case: false,
        canonical: true,
        sketch_size: sketch_params.expected_size() as u32,
        hash_type: String::from(hash_type),
        hash_bits,
        hash_seed,
        sketches: sketches?,
    })
}

pub fn sketch_stream<'a>(
    reader: Box<dyn Read + 'a>,
    name: &str,
    sketch_params: &SketchParams,
    filters: &FilterParams,
) -> Result<Sketch> {
    let mut filters = filters.clone();
    let mut sketcher = sketch_params.create_sketcher();
    parse_sequence_reader(
        reader,
        |seq_type| {
            // disable filtering for FASTA files unless it was explicitly specified
            if filters.filter_on.is_none() {
                filters.filter_on = match seq_type {
                    "FASTA" => Some(false),
                    "FASTQ" => Some(true),
                    _ => panic!("Unknown sequence type"),
                };
            }
        },
        |seq| {
            sketcher.process(seq);
        },
    )
    .map_err(|e| format_err!("{}", e.to_string()))?;

    let (seq_len, n_kmers) = sketcher.total_bases_and_kmers();
    let hashes = sketcher.to_vec();

    // do filtering
    let (mut filtered_hashes, low_abun) = filters.filter_sketch(&hashes);
    let mut filter_stats = filters.serialize_filter_params();
    if let Some(la) = low_abun {
        filter_stats.insert(String::from("minCopies"), la.to_string());
    }

    sketch_params.process_post_filter(&mut filtered_hashes, name)?;

    Ok(Sketch::new(
        name,
        seq_len,
        n_kmers,
        filtered_hashes,
        &filter_stats,
    ))
}
