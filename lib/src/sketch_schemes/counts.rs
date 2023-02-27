use needletail::bitkmer::{bitmer_to_bytes, reverse_complement};
use needletail::Sequence;

use crate::sketch_schemes::{KmerCount, SketchParams, SketchScheme};

#[derive(Clone)]
pub struct AllCountsSketcher {
    counts: Vec<u32>,
    total_bases: u64,
    k: u8,
}

impl AllCountsSketcher {
    pub fn new(k: u8) -> Self {
        // TODO: should we take a size parameter or the like and clip this?
        AllCountsSketcher {
            counts: vec![0; 4usize.pow(k.into())],
            total_bases: 0,
            k,
        }
    }
}

impl SketchScheme for AllCountsSketcher {
    fn process<'s, 'a: 's, 'b>(&'a mut self, seq: &'s dyn Sequence<'b>)
    where
        's: 'b,
    {
        for (_, kmer, _) in seq.normalize(false).bit_kmers(self.k, false) {
            self.counts[kmer.0 as usize] = self.counts[kmer.0 as usize].saturating_add(1);
        }
    }

    fn total_bases_and_kmers(&self) -> (u64, u64) {
        (
            self.total_bases,
            self.counts.iter().map(|x| u64::from(*x)).sum(),
        )
    }

    fn to_vec(&self) -> Vec<KmerCount> {
        let mut counts = self.counts.clone();
        let mut results = Vec::with_capacity(self.counts.len());
        for ix in 0u64..counts.len() as u64 {
            let mut count = counts[ix as usize];
            if count == 0 {
                continue;
            }
            let extra_count = self.counts[reverse_complement((ix, self.k)).0 as usize];
            counts[reverse_complement((ix, self.k)).0 as usize] = 0;
            count += extra_count;
            let new_item = KmerCount {
                hash: ix,
                kmer: bitmer_to_bytes((ix, self.k)),
                count,
                extra_count,
                label: None,
            };
            results.push(new_item);
        }
        results
    }

    fn parameters(&self) -> SketchParams {
        SketchParams::AllCounts {
            kmer_length: self.k,
        }
    }
}
