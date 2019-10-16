// Compile the capnp schema with `capnp compile ./src/serialization/finch.capnp -orust`
//
// paths are broken so you need `%s/::finch_capnp::/super::/g`; see:
// https://github.com/capnproto/capnproto-rust/issues/16
//
// also add the following to disable most warnings:
// #![allow(clippy::all)]
// #![allow(dead_code)]
mod finch_capnp;

mod json;
mod mash;
mod mash_capnp;

use std::collections::HashMap;
use std::io::{BufRead, Write};

use capnp::message;
use capnp::serialize as capnp_serialize;

use crate::filtering::FilterParams;
use crate::serialization::finch_capnp::{multisketch, sketch_params, SketchMethod};
pub use crate::serialization::json::{MultiSketch, Sketch};
pub use crate::serialization::mash::{read_mash_file, write_mash_file};
use crate::sketch_schemes::{KmerCount, SketchParams};
use crate::Result;

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

// TODO: use this as the "base" sketch model within finch some day
#[derive(Clone, Debug)]
pub struct FinchSketch {
    pub name: String,
    pub seq_length: u64,
    pub num_valid_kmers: u64,
    pub comment: String,

    pub hashes: Vec<KmerCount>,
    pub filter_params: FilterParams,
    pub sketch_params: SketchParams,
}

impl FinchSketch {
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }
}

impl Into<Sketch> for FinchSketch {
    fn into(self) -> Sketch {
        Sketch::new(
            &self.name,
            self.seq_length,
            self.num_valid_kmers,
            self.hashes,
            &self.filter_params.to_serialized(),
        )
    }
}

impl Into<Vec<FinchSketch>> for &MultiSketch {
    fn into(self) -> Vec<FinchSketch> {
        let empty_hashmap = HashMap::new();
        let mut sketches = Vec::with_capacity(self.sketches.len());
        let sketch_params = self.get_params().unwrap();
        for sketch in &self.sketches {
            let filters = sketch.filters.as_ref().unwrap_or(&empty_hashmap);
            let filter_params = FilterParams::from_serialized(filters);
            sketches.push(FinchSketch {
                name: sketch.name.clone(),
                seq_length: sketch.seq_length.unwrap_or(0),
                num_valid_kmers: sketch.num_valid_kmers.unwrap_or(0),
                comment: sketch.comment.clone().unwrap_or_else(|| "".to_string()),
                hashes: sketch.hashes.clone(),
                filter_params,
                sketch_params: sketch_params.clone(),
            });
        }
        sketches
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

fn get_sketch_params(cap_sketch_params: sketch_params::Reader) -> Result<SketchParams> {
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

fn saturated_u64_to_u32(v: u64) -> u32 {
    if v > u64::from(std::u32::MAX) {
        std::u32::MAX
    } else {
        v as u32
    }
}

pub fn write_finch_file(mut file: &mut dyn Write, sketches: &[FinchSketch]) -> Result<()> {
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
            cap_hash.set_count(saturated_u64_to_u32(hash.count));
            cap_hash.set_extra_count(saturated_u64_to_u32(hash.extra_count));
        }

        let mut cap_filter_params = cap_sketch.reborrow().init_filter_params();
        cap_filter_params.set_filtered(sketch.filter_params.filter_on.unwrap_or(false));
        cap_filter_params.set_low_abun_filter(sketch.filter_params.abun_filter.0.unwrap_or(0));
        cap_filter_params.set_high_abun_filter(
            sketch
                .filter_params
                .abun_filter
                .1
                .unwrap_or(::std::u64::MAX),
        );
        cap_filter_params.set_err_filter(sketch.filter_params.err_filter);
        cap_filter_params.set_strand_filter(sketch.filter_params.strand_filter);

        let sketch_params = &sketch.sketch_params;
        let cap_sketch_params = cap_sketch.reborrow().init_sketch_params();
        set_sketch_params(cap_sketch_params, &sketch_params);
    }

    capnp_serialize::write_message(&mut file, &message)?;
    Ok(())
}

pub fn read_finch_file(mut file: &mut dyn BufRead) -> Result<MultiSketch> {
    let options = *message::ReaderOptions::new().traversal_limit_in_words(64 * 1024 * 1024);
    let reader = capnp_serialize::read_message(&mut file, options)?;
    let cap_data: multisketch::Reader = reader.get_root::<multisketch::Reader>()?;
    let cap_sketches = cap_data.get_sketches()?;

    let mut sketches = Vec::with_capacity(cap_sketches.len() as usize);
    for cap_sketch in cap_sketches {
        let cap_hashes = cap_sketch.get_hashes()?;
        let mut hashes = Vec::with_capacity(cap_hashes.len() as usize);
        for cap_hash in cap_hashes {
            hashes.push(KmerCount {
                hash: cap_hash.get_hash(),
                kmer: cap_hash.get_kmer()?.to_vec(),
                count: u64::from(cap_hash.get_count()),
                extra_count: u64::from(cap_hash.get_extra_count()),
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
            ::std::u64::MAX => None,
            i => Some(i),
        };
        let filter_params = FilterParams {
            filter_on: Some(cap_filter_params.get_filtered()),
            abun_filter: (low_abun_filter, high_abun_filter),
            err_filter: cap_filter_params.get_err_filter(),
            strand_filter: cap_filter_params.get_strand_filter(),
        };

        sketches.push(FinchSketch {
            name: cap_sketch.get_name()?.to_string(),
            seq_length: cap_sketch.get_seq_length(),
            num_valid_kmers: cap_sketch.get_num_valid_kmers(),
            comment: cap_sketch.get_comment()?.to_string(),

            hashes,
            sketch_params,
            filter_params,
        });
    }

    let json_sketches: Vec<Sketch> = sketches.iter().map(|x| (*x).clone().into()).collect();
    let sketch_params = &sketches[0].sketch_params;
    let (hash_type, hash_bits, hash_seed, scale) = sketch_params.hash_info();
    Ok(MultiSketch {
        alphabet: "ACGT".to_string(),
        preserve_case: false,
        canonical: true,

        sketch_size: sketch_params.expected_size() as u32,
        kmer: sketch_params.k(),
        hash_type: hash_type.to_string(),
        hash_bits,
        hash_seed,
        scale,
        sketches: json_sketches,
    })
}
