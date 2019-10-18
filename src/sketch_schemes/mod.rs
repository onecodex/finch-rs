pub mod counts;
mod hashing;
pub mod mash;
pub mod scaled;

use failure::bail;
use needletail::SequenceRecord;

use crate::Result;
pub use hashing::ItemHash;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Hash, Serialize)]
pub struct KmerCount {
    pub hash: ItemHash,
    pub kmer: Vec<u8>,
    pub count: u32,
    pub extra_count: u32,
}

pub trait SketchScheme {
    fn process(&mut self, seq: SequenceRecord);
    fn total_bases_and_kmers(&self) -> (u64, u64);
    fn to_vec(&self) -> Vec<KmerCount>;
}

#[derive(Clone, Debug)]
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

    pub fn process_post_filter(&self, kmers: &mut Vec<KmerCount>, name: &str) -> Result<()> {
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

    // TODO: check compatibility with another param set?
}
