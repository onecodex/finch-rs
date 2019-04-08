#![cfg_attr(feature = "python", feature(specialization))]

#[cfg(feature = "mash_format")]
extern crate capnp;
#[macro_use] extern crate serde_derive;

use std::io::{Read, Seek};
use std::path::Path;
use std::result::Result as StdResult;
use failure::{Error, bail, format_err};
use needletail::fastx::{fastx_cli, fastx_stream};

use crate::filtering::{FilterParams, filter_sketch};
use crate::minhashes::MinHashKmers;
use crate::serialization::{Sketch, MultiSketch};

pub mod minhashes;
pub mod filtering;
pub mod distance;
pub mod serialization;
pub mod statistics;
#[cfg(feature = "mash_format")]
mod mash_capnp;
#[cfg(feature = "python")]
pub mod python;

pub type Result<T> = StdResult<T, Error>;

pub fn mash_files(filenames: &[&str], n_hashes: usize, final_size: usize, kmer_length: u8, filters: &mut FilterParams, no_strict: bool, seed: u64) -> Result<MultiSketch> {
    let mut sketches = Vec::with_capacity(filenames.len());
    for filename in filenames {
        let mut seq_len = 0u64;
        let path = Path::new(filename);
        let mut minhash = match filters.filter_on {
            Some(true) | None => MinHashKmers::new(n_hashes, seed),
            Some(false) => MinHashKmers::new(final_size, seed),
        };
        fastx_cli(
            path.to_str()
                .ok_or(format_err!("Couldn't make path into string"))?,
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
                for (_, kmer, is_rev_complement) in seq.normalize(false).kmers(kmer_length, true) {
                    let rc_count = if is_rev_complement { 1u8 } else { 0u8 };
                    minhash.push(kmer, rc_count);
                }
            },
        )
        .map_err(|e| format_err!("{}", e.to_string()))?;

        let n_kmers = minhash.total_kmers() as u64;
        let hashes = minhash.into_vec();
        let (mut filtered_hashes, filter_stats) = filter_sketch(&hashes, &filters);
        filtered_hashes.truncate(final_size);
        if !no_strict && filtered_hashes.len() < final_size {
            bail!("{} had too few kmers ({}) to sketch", filename, filtered_hashes.len());
        }

        // directory should be clipped from filename
        let basename = path.file_name().ok_or(format_err!("Couldn't get filename from path"))?;
        let sketch = Sketch::new(basename.to_str().ok_or(format_err!("Couldn't make filename into string"))?,
                                     seq_len, n_kmers, filtered_hashes, &filter_stats);
        sketches.push(sketch);
    }
    Ok(MultiSketch {
        kmer: kmer_length,
        alphabet: String::from("ACGT"),
        preserveCase: false,
        canonical: true,
        sketchSize: final_size as u32,
        hashType: String::from("MurmurHash3_x64_128"),
        hashBits: 64u16,
        hashSeed: seed,
        sketches,
    })
}


pub fn mash_stream<R>(reader: R, n_hashes: usize, final_size: usize, kmer_length: u8,
                      filters: &mut FilterParams, no_strict: bool, seed: u64) -> Result<Sketch> where
    R: Read + Seek,
{
    let mut seq_len = 0u64;
    let mut minhash = match filters.filter_on {
        Some(true) | None => MinHashKmers::new(n_hashes, seed),
        Some(false) => MinHashKmers::new(final_size, seed),
    };
    fastx_stream(
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
        }, |seq| {
            seq_len += seq.seq.len() as u64;
            for (_, kmer, is_rev_complement) in seq.normalize(false).kmers(kmer_length, true) {
                let rc_count = if is_rev_complement { 1u8 } else { 0u8 };
                minhash.push(kmer, rc_count);
            }
        }).map_err(|e| format_err!("{}", e.to_string()))?;

    let n_kmers = minhash.total_kmers() as u64;
    let hashes = minhash.into_vec();
    let (mut filtered_hashes, filter_stats) = filter_sketch(&hashes, &filters);
    filtered_hashes.truncate(final_size);
    if !no_strict && filtered_hashes.len() < final_size {
        bail!(
            "Stream had too few kmers ({}) to sketch",
            filtered_hashes.len()
        );
    }

        Ok(Sketch::new("", seq_len, n_kmers, filtered_hashes, &filter_stats))
}
