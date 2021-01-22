pub mod counts;
mod hashing;
pub mod mash;
pub mod scaled;

use needletail::parser::SequenceRecord;
use serde::{Deserialize, Serialize};

use crate::bail;
use crate::errors::FinchResult;
use crate::filtering::FilterParams;
use crate::serialization::Sketch;
pub use hashing::ItemHash;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Hash, Serialize)]
pub struct KmerCount {
    pub hash: ItemHash,
    pub kmer: Vec<u8>,
    pub count: u32,
    pub extra_count: u32,
    pub label: Option<Vec<u8>>,
}

pub trait SketchScheme {
    fn process(&mut self, seq: SequenceRecord);
    fn total_bases_and_kmers(&self) -> (u64, u64);
    fn to_vec(&self) -> Vec<KmerCount>;
    fn parameters(&self) -> SketchParams;

    fn to_sketch(&self) -> Sketch {
        // TODO: maybe this should be the primary teardown method for
        // sketching and sketch_stream should wrap it?
        // TODO: this doesn't really use filtering
        // TODO: also the pass-through for the post-filtering trimming is
        // weird for SketchParams::Mash
        let (seq_length, num_valid_kmers) = self.total_bases_and_kmers();
        let hashes = self.to_vec();
        Sketch {
            name: "".to_string(),
            seq_length,
            num_valid_kmers,
            comment: "".to_string(),
            hashes,
            filter_params: FilterParams::default(),
            sketch_params: self.parameters(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum SketchParams {
    Mash {
        kmers_to_sketch: usize,
        final_size: usize,
        no_strict: bool,
        kmer_length: u8,
        hash_seed: u64,
    },
    Scaled {
        kmers_to_sketch: usize,
        kmer_length: u8,
        scale: f64,
        hash_seed: u64,
    },
    AllCounts {
        kmer_length: u8,
    },
}

impl Default for SketchParams {
    fn default() -> Self {
        SketchParams::Mash {
            kmers_to_sketch: 1000,
            final_size: 1000,
            no_strict: false,
            kmer_length: 21,
            hash_seed: 0,
        }
    }
}

impl SketchParams {
    pub fn create_sketcher(&self) -> Box<dyn SketchScheme> {
        match self {
            SketchParams::Mash {
                kmers_to_sketch,
                kmer_length,
                hash_seed,
                ..
            } => Box::new(mash::MashSketcher::new(
                *kmers_to_sketch,
                *kmer_length,
                *hash_seed,
            )),
            SketchParams::Scaled {
                kmers_to_sketch,
                kmer_length,
                scale,
                hash_seed,
            } => Box::new(scaled::ScaledSketcher::new(
                *kmers_to_sketch,
                *scale,
                *kmer_length,
                *hash_seed,
            )),
            SketchParams::AllCounts { kmer_length } => {
                Box::new(counts::AllCountsSketcher::new(*kmer_length))
            }
        }
    }

    pub fn process_post_filter(&self, kmers: &mut Vec<KmerCount>, name: &str) -> FinchResult<()> {
        if let SketchParams::Mash {
            final_size,
            no_strict,
            ..
        } = self
        {
            kmers.truncate(*final_size);
            if !no_strict && kmers.len() < *final_size {
                bail!("{} had too few kmers ({}) to sketch", name, kmers.len(),);
            }
        }
        Ok(())
    }

    pub fn k(&self) -> u8 {
        match self {
            SketchParams::Mash { kmer_length, .. } => *kmer_length,
            SketchParams::Scaled { kmer_length, .. } => *kmer_length,
            SketchParams::AllCounts { kmer_length, .. } => *kmer_length,
        }
    }

    pub fn hash_info(&self) -> (&str, u16, u64, Option<f64>) {
        match self {
            SketchParams::Mash { hash_seed, .. } => ("MurmurHash3_x64_128", 64, *hash_seed, None),
            SketchParams::Scaled {
                hash_seed, scale, ..
            } => ("MurmurHash3_x64_128", 64, *hash_seed, Some(*scale)),
            SketchParams::AllCounts { .. } => ("None", 0, 0, None),
        }
    }

    pub fn expected_size(&self) -> usize {
        match self {
            SketchParams::Mash { final_size, .. } => *final_size,
            SketchParams::Scaled {
                kmers_to_sketch, ..
            } => *kmers_to_sketch,
            SketchParams::AllCounts { kmer_length, .. } => 4usize.pow(u32::from(*kmer_length)),
        }
    }

    pub fn from_sketches(sketches: &[Sketch]) -> FinchResult<Self> {
        let first_params = sketches[0].sketch_params.clone();
        for (ix, sketch) in sketches.iter().enumerate().skip(1) {
            let params = &sketch.sketch_params;
            if let Some((mismatched_param, v1, v2)) = first_params.check_compatibility(&params) {
                bail!(
                    "First sketch has {} {}, but sketch {} has {0} {}",
                    mismatched_param,
                    v1,
                    ix + 1,
                    v2,
                );
            }
            // TODO: harmonize scaled/non-scaled sketches?
            // TODO: harminize sketch sizes?
            // TODO: do something with no_strict and final_size
        }

        Ok(first_params)
    }

    /// Return any sketch parameter difference that would make comparisons
    /// between sketches generated by these parameter sets not work.
    ///
    /// Note this doesn't actually check the enum variants themselves, but it
    /// should still break if there are different variants because the hash
    /// types should be different.
    pub fn check_compatibility(&self, other: &SketchParams) -> Option<(&str, String, String)> {
        if self.k() != other.k() {
            return Some(("k", self.k().to_string(), other.k().to_string()));
        }
        if self.hash_info().0 != other.hash_info().0 {
            return Some((
                "hash type",
                self.hash_info().0.to_string(),
                other.hash_info().0.to_string(),
            ));
        }
        if self.hash_info().1 != other.hash_info().1 {
            return Some((
                "hash bits",
                self.hash_info().1.to_string(),
                other.hash_info().1.to_string(),
            ));
        }
        if self.hash_info().2 != other.hash_info().2 {
            return Some((
                "hash seed",
                self.hash_info().2.to_string(),
                other.hash_info().2.to_string(),
            ));
        }

        None
    }
}
