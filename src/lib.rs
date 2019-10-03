#![cfg_attr(feature = "python", feature(specialization))]

#[cfg(feature = "mash_format")]
extern crate capnp;
#[macro_use]
extern crate serde_derive;

use std::fs::File;
use std::io::{stdin, Read};
use std::path::Path;
use std::result::Result as StdResult;

use failure::{bail, format_err, Error};
use needletail::formats::parse_sequence_reader;
use rayon::prelude::*;

use crate::filtering::FilterParams;
use crate::hash_schemes::minhashes::MinHashKmers;
use crate::hash_schemes::HashScheme;
use crate::serialization::{MultiSketch, Sketch};

pub mod distance;
pub mod filtering;
pub mod hash_schemes;
// TODO: it would be nice if there was a `pub(in main)` or something for main_parsing?
pub mod main_parsing;
#[cfg(feature = "mash_format")]
mod mash_capnp;
#[cfg(feature = "python")]
pub mod python;
pub mod serialization;
pub mod statistics;

pub type Result<T> = StdResult<T, Error>;

pub fn mash_files(
    filenames: &[&str],
    n_hashes: usize,
    final_size: usize,
    kmer_length: u8,
    filters: &FilterParams,
    no_strict: bool,
    seed: u64,
) -> Result<MultiSketch> {
    let sketches: Result<Vec<Sketch>> = filenames
        .par_iter()
        .map(|filename| {
            let sketch = if filename == &"-" {
                // special case for stdin
                let sin = stdin();
                mash_stream(
                    sin.lock(),
                    filename,
                    n_hashes,
                    final_size,
                    kmer_length,
                    &filters,
                    no_strict,
                    seed,
                )?
            } else {
                mash_stream(
                    File::open(&Path::new(filename))?,
                    filename,
                    n_hashes,
                    final_size,
                    kmer_length,
                    &filters,
                    no_strict,
                    seed,
                )?
            };
            Ok(sketch)
        })
        .collect();
    Ok(MultiSketch {
        kmer: kmer_length,
        alphabet: String::from("ACGT"),
        preserveCase: false,
        canonical: true,
        sketchSize: final_size as u32,
        hashType: String::from("MurmurHash3_x64_128"),
        hashBits: 64u16,
        hashSeed: seed,
        sketches: sketches?,
    })
}

pub fn mash_stream<R: Read>(
    reader: R,
    name: &str,
    n_hashes: usize,
    final_size: usize,
    kmer_length: u8,
    filters: &FilterParams,
    no_strict: bool,
    seed: u64,
) -> Result<Sketch> {
    let mut seq_len = 0u64;
    let mut filters = filters.clone();
    let mut minhash = match filters.filter_on {
        Some(true) | None => MinHashKmers::new(n_hashes, kmer_length, seed),
        Some(false) => MinHashKmers::new(final_size, kmer_length, seed),
    };
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
            seq_len += seq.seq.len() as u64;
            minhash.process(seq);
        },
    )
    .map_err(|e| format_err!("{}", e.to_string()))?;

    let n_kmers = minhash.total_kmers() as u64;
    let hashes = minhash.into_vec();

    // do filtering
    let (mut filtered_hashes, low_abun) = filters.filter_sketch(&hashes);
    let mut filter_stats = filters.serialize_filter_params();
    if let Some(la) = low_abun {
        filter_stats.insert(String::from("minCopies"), la.to_string());
    }

    filtered_hashes.truncate(final_size);
    if !no_strict && filtered_hashes.len() < final_size {
        bail!(
            "{} had too few kmers ({}) to sketch",
            name,
            filtered_hashes.len(),
        );
    }

    Ok(Sketch::new(
        name,
        seq_len,
        n_kmers,
        filtered_hashes,
        &filter_stats,
    ))
}
