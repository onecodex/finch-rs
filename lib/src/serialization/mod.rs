// Those 2 files are generated so we just ignore everything in terms
#[allow(clippy::all)]
#[allow(dead_code)]
#[cfg_attr(rustfmt, rustfmt_skip)]
mod finch_capnp;
#[allow(clippy::all)]
#[allow(dead_code)]
#[cfg_attr(rustfmt, rustfmt_skip)]
mod mash_capnp;

mod json;
mod mash;

use std::io::{BufRead, Write};

use capnp::message;
use capnp::serialize as capnp_serialize;
use serde::{Deserialize, Serialize};

use crate::errors::FinchResult;
use crate::filtering::FilterParams;
use crate::serialization::finch_capnp::{multisketch, sketch_params, SketchMethod};
pub use crate::serialization::json::{JsonSketch, MultiSketch};
pub use crate::serialization::mash::{read_mash_file, write_mash_file};
use crate::sketch_schemes::{KmerCount, SketchParams};

pub const FINCH_EXT: &str = ".sk";
pub const FINCH_BIN_EXT: &str = ".bsk";
pub const MASH_EXT: &str = ".msh";

#[derive(Debug, Serialize, Deserialize)]
pub struct SketchDistance {
    pub containment: f64,
    pub jaccard: f64,
    #[serde(rename = "mashDistance")]
    pub mash_distance: f64,
    #[serde(rename = "commonHashes")]
    pub common_hashes: u64,
    #[serde(rename = "totalHashes")]
    pub total_hashes: u64,
    pub query: String,
    pub reference: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Sketch {
    pub name: String,
    pub seq_length: u64,
    pub num_valid_kmers: u64,
    pub comment: String,

    pub hashes: Vec<KmerCount>,
    pub filter_params: FilterParams,
    pub sketch_params: SketchParams,
}

impl Sketch {
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }
}

fn set_sketch_params(mut cap_sketch_params: sketch_params::Builder, sketch_params: &SketchParams) {
    match *sketch_params {
        SketchParams::Mash {
            kmers_to_sketch,
            final_size,
            no_strict,
            kmer_length,
            hash_seed,
        } => {
            cap_sketch_params.set_sketch_method(SketchMethod::MurmurHash3);
            cap_sketch_params.set_kmer_length(kmer_length);
            cap_sketch_params.set_kmers_to_sketch(kmers_to_sketch as u64);
            cap_sketch_params.set_hash_seed(hash_seed);
            cap_sketch_params.set_final_size(final_size as u64);
            cap_sketch_params.set_no_strict(no_strict);
        }
        SketchParams::Scaled {
            kmers_to_sketch,
            kmer_length,
            scale,
            hash_seed,
        } => {
            cap_sketch_params.set_sketch_method(SketchMethod::MurmurHash3Scaled);
            cap_sketch_params.set_kmer_length(kmer_length);
            cap_sketch_params.set_kmers_to_sketch(kmers_to_sketch as u64);
            cap_sketch_params.set_hash_seed(hash_seed);
            cap_sketch_params.set_scale(scale);
        }
        SketchParams::AllCounts { kmer_length } => {
            cap_sketch_params.set_sketch_method(SketchMethod::None);
            cap_sketch_params.set_kmer_length(kmer_length);
        }
    }
}

fn get_sketch_params(cap_sketch_params: sketch_params::Reader) -> FinchResult<SketchParams> {
    Ok(match cap_sketch_params.get_sketch_method()? {
        SketchMethod::MurmurHash3 => SketchParams::Mash {
            kmers_to_sketch: cap_sketch_params.get_kmers_to_sketch() as usize,
            final_size: cap_sketch_params.get_final_size() as usize,
            no_strict: cap_sketch_params.get_no_strict(),
            kmer_length: cap_sketch_params.get_kmer_length(),
            hash_seed: cap_sketch_params.get_hash_seed(),
        },
        SketchMethod::MurmurHash3Scaled => SketchParams::Scaled {
            kmers_to_sketch: cap_sketch_params.get_kmers_to_sketch() as usize,
            kmer_length: cap_sketch_params.get_kmer_length(),
            scale: cap_sketch_params.get_scale(),
            hash_seed: cap_sketch_params.get_hash_seed(),
        },
        SketchMethod::None => SketchParams::AllCounts {
            kmer_length: cap_sketch_params.get_kmer_length(),
        },
    })
}

pub fn write_finch_file(file: &mut dyn Write, sketches: &[Sketch]) -> FinchResult<()> {
    let mut message = message::Builder::new_default();
    let finch_file: multisketch::Builder = message.init_root::<multisketch::Builder>();

    let mut cap_sketches = finch_file.init_sketches(sketches.len() as u32);
    for (i, sketch) in sketches.iter().enumerate() {
        let mut cap_sketch = cap_sketches.reborrow().get(i as u32);
        cap_sketch.set_name(&sketch.name);
        cap_sketch.set_seq_length(sketch.seq_length);
        cap_sketch.set_num_valid_kmers(sketch.num_valid_kmers);
        cap_sketch.set_comment(&sketch.comment);

        // TODO: we should probably error if hashes.len() > 2**32?
        // (and handle these `as u32`s a little better in general
        let mut hashes = cap_sketch
            .reborrow()
            .init_hashes(sketch.hashes.len() as u32);
        for (j, hash) in sketch.hashes.iter().enumerate() {
            let mut cap_hash = hashes.reborrow().get(j as u32);
            cap_hash.set_hash(hash.hash);
            cap_hash.set_kmer(&hash.kmer);
            cap_hash.set_count(hash.count);
            cap_hash.set_extra_count(hash.extra_count);
            if let Some(label) = &hash.label {
                cap_hash.set_label(label);
            }
        }

        let mut cap_filter_params = cap_sketch.reborrow().init_filter_params();
        cap_filter_params.set_filtered(sketch.filter_params.filter_on.unwrap_or(false));
        cap_filter_params.set_low_abun_filter(sketch.filter_params.abun_filter.0.unwrap_or(0));
        cap_filter_params
            .set_high_abun_filter(sketch.filter_params.abun_filter.1.unwrap_or(u32::MAX));
        cap_filter_params.set_err_filter(sketch.filter_params.err_filter);
        cap_filter_params.set_strand_filter(sketch.filter_params.strand_filter);

        let sketch_params = &sketch.sketch_params;
        let cap_sketch_params = cap_sketch.reborrow().init_sketch_params();
        set_sketch_params(cap_sketch_params, sketch_params);
    }

    capnp_serialize::write_message(file, &message)?;
    Ok(())
}

pub fn read_finch_file(file: &mut dyn BufRead) -> FinchResult<Vec<Sketch>> {
    let options = *message::ReaderOptions::new().traversal_limit_in_words(Some(1024 * 1024 * 1024));
    let reader = capnp_serialize::read_message(file, options)?;
    let cap_data: multisketch::Reader = reader.get_root::<multisketch::Reader>()?;
    let cap_sketches = cap_data.get_sketches()?;

    let mut sketches = Vec::with_capacity(cap_sketches.len() as usize);
    for cap_sketch in cap_sketches {
        let cap_hashes = cap_sketch.get_hashes()?;
        let mut hashes = Vec::with_capacity(cap_hashes.len() as usize);
        for cap_hash in cap_hashes {
            let label = if cap_hash.has_label() {
                Some(cap_hash.get_label()?.to_vec())
            } else {
                None
            };
            hashes.push(KmerCount {
                hash: cap_hash.get_hash(),
                kmer: cap_hash.get_kmer()?.to_vec(),
                count: cap_hash.get_count(),
                extra_count: cap_hash.get_extra_count(),
                label,
            });
        }

        let cap_sketch_params = cap_sketch.get_sketch_params()?;
        let sketch_params = get_sketch_params(cap_sketch_params)?;

        let cap_filter_params = cap_sketch.get_filter_params()?;
        let low_abun_filter = match cap_filter_params.get_low_abun_filter() {
            0 => None,
            i => Some(i),
        };
        let high_abun_filter = match cap_filter_params.get_high_abun_filter() {
            ::std::u32::MAX => None,
            i => Some(i),
        };
        let filter_params = FilterParams {
            filter_on: Some(cap_filter_params.get_filtered()),
            abun_filter: (low_abun_filter, high_abun_filter),
            err_filter: cap_filter_params.get_err_filter(),
            strand_filter: cap_filter_params.get_strand_filter(),
        };

        sketches.push(Sketch {
            name: cap_sketch.get_name()?.to_string(),
            seq_length: cap_sketch.get_seq_length(),
            num_valid_kmers: cap_sketch.get_num_valid_kmers(),
            comment: cap_sketch.get_comment()?.to_string(),

            hashes,
            sketch_params,
            filter_params,
        });
    }
    Ok(sketches)
}
