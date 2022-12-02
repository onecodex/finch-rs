use std::fs::File;
use std::io::{stdin, BufReader, Read};
use std::path::Path;

use memmap::MmapOptions;
use needletail::parse_fastx_reader;
use needletail::parser::Format;
use rayon::prelude::*;

use crate::filtering::FilterParams;
use crate::serialization::{
    read_finch_file, read_mash_file, MultiSketch, Sketch, FINCH_BIN_EXT, FINCH_EXT, MASH_EXT,
};
use crate::sketch_schemes::SketchParams;

pub mod distance;
pub mod filtering;
pub mod sketch_schemes;
// it would be nice if there was a `pub(in main)` or something for
// main_parsing so we don't import it for `lib` itself
pub mod errors;
#[cfg(feature = "python")]
pub mod python;
pub mod serialization;
pub mod statistics;

use crate::errors::FinchResult;

pub fn sketch_files(
    filenames: &[&str],
    sketch_params: &SketchParams,
    filters: &FilterParams,
) -> FinchResult<Vec<Sketch>> {
    let sketches: FinchResult<Vec<Sketch>> = filenames
        .par_iter()
        .map(|filename| {
            // open the file with a special case to handle stdin
            let reader: Box<dyn Read + Send> = if filename == &"-" {
                // We're not locking it so technically not thread safe
                Box::new(stdin())
            } else {
                Box::new(File::open(&Path::new(filename))?)
            };
            // sketch!
            sketch_stream(reader, filename, sketch_params, filters)
        })
        .collect();
    sketches
}

pub fn sketch_stream<'a>(
    reader: Box<dyn Read + Send + 'a>,
    name: &str,
    sketch_params: &SketchParams,
    filters: &FilterParams,
) -> FinchResult<Sketch> {
    let mut filter_params = filters.clone();
    let mut sketcher = sketch_params.create_sketcher();
    // TODO: remove expects after removing failure
    let mut fastx_reader = parse_fastx_reader(reader).expect("valid file TODO");
    let mut seq_type = None;
    while let Some(record) = fastx_reader.next() {
        let seqrec = record.expect("invalid record");
        if seq_type.is_none() {
            seq_type = Some(seqrec.format());
        }
        sketcher.process(&seqrec);
    }

    // disable filtering for FASTA files unless it was explicitly specified
    if filter_params.filter_on.is_none() {
        filter_params.filter_on = match seq_type.expect("Should have got a type") {
            Format::Fasta => Some(false),
            Format::Fastq => Some(true),
        };
    }

    let (seq_length, num_valid_kmers) = sketcher.total_bases_and_kmers();
    let hashes = sketcher.to_vec();

    // do filtering
    let mut filtered_hashes = filter_params.filter_counts(&hashes);
    sketch_params.process_post_filter(&mut filtered_hashes, name)?;

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

pub fn open_sketch_file<P: AsRef<Path>>(path: P) -> FinchResult<Vec<Sketch>> {
    let p = path.as_ref();
    let filename = p
        .file_name()
        .ok_or_else(|| format_err!("Path does not have a filename: {:?}", p))?
        .to_string_lossy();
    let file = File::open(p).map_err(|_| format_err!("Error opening {:?}", p))?;
    if filename.ends_with(MASH_EXT) {
        let mut buf_reader = BufReader::new(file);
        read_mash_file(&mut buf_reader)
    } else if filename.ends_with(FINCH_BIN_EXT) {
        let mut buf_reader = BufReader::new(file);
        read_finch_file(&mut buf_reader)
    } else if filename.ends_with(FINCH_EXT) || filename.ends_with(".json") {
        let mapped = unsafe { MmapOptions::new().map(&file)? };
        let multisketch: MultiSketch =
            serde_json::from_slice(&mapped).map_err(|_| format_err!("Error parsing {:?}", &p))?;
        multisketch.to_sketches()
    } else {
        Err(format_err!("File suffix is not *.bsk, *.msh, or *.sk"))
    }
}
