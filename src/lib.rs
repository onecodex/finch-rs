extern crate needletail;
extern crate murmurhash3;
extern crate serde;
#[macro_use] extern crate serde_derive;
extern crate serde_json;

use std::path::Path;
use needletail::fastx::fastx_cli;

use filtering::{FilterParams, filter_sketch};
use minhashes::MinHashKmers;
use serialization::{JSONSketch, JSONMultiSketch};

pub mod minhashes;
pub mod filtering;
pub mod distance;
pub mod serialization;
pub mod statistics;

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
