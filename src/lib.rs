#![cfg_attr(feature = "python", feature(specialization))]

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
    Ok(MultiSketch::from_sketches(&sketches?))
}

pub fn sketch_stream<'a>(
    reader: Box<dyn Read + 'a>,
    name: &str,
    sketch_params: &SketchParams,
    filters: &FilterParams,
) -> Result<Sketch> {
    let mut filter_params = filters.clone();
    let mut sketcher = sketch_params.create_sketcher();
    parse_sequence_reader(
        reader,
        |seq_type| {
            // disable filtering for FASTA files unless it was explicitly specified
            if filter_params.filter_on.is_none() {
                filter_params.filter_on = match seq_type {
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

    let (seq_length, num_valid_kmers) = sketcher.total_bases_and_kmers();
    let hashes = sketcher.to_vec();

    // do filtering
    let mut filtered_hashes = filter_params.filter_sketch(&hashes);
    sketch_params.process_post_filter(&mut filtered_hashes, name)?;
    // let filter_stats = filters.to_serialized();

    Ok(Sketch {
        name: name.to_string(),
        seq_length,
        num_valid_kmers,
        comment: "".to_string(),
        hashes: filtered_hashes,
        filter_params,
        sketch_params: sketch_params.clone(),
    })
}
