use needletail::bitkmer::{bitmer_to_bytes, reverse_complement};
use needletail::{Sequence, SequenceRecord};

use crate::hash_schemes::{HashScheme, KmerCount};

#[derive(Clone)]
pub struct AllKmers {
    counts: Vec<u16>,
    k: u8,
}

impl AllKmers {
    pub fn new(k: u8) -> Self {
        // TODO: should we take a size parameter or the like and clip this?
        AllKmers {
            counts: vec![0; 4usize.pow(k.into())],
            k,
        }
    }
}

impl HashScheme for AllKmers {
    fn process(&mut self, seq: SequenceRecord) {
        for (_, kmer, _) in seq.normalize(false).bit_kmers(self.k, false) {
            self.counts[kmer.0 as usize] += 1;
        }
    }

    fn total_kmers(&self) -> usize {
        let total: u16 = self.counts.iter().sum();
        total as usize
    }

    fn into_vec(mut self) -> Vec<KmerCount> {
        let mut results = Vec::with_capacity(self.counts.len());
        for ix in 0u64..self.counts.len() as u64 {
            let mut count = self.counts[ix as usize];
            if count == 0 {
                continue;
            }
            let extra_count = self.counts[reverse_complement((ix, self.k)).0 as usize];
            self.counts[reverse_complement((ix, self.k)).0 as usize] = 0;
            count += extra_count;
            let new_item = KmerCount {
                hash: ix,
                kmer: bitmer_to_bytes((ix as u64, self.k)),
                count,
                extra_count,
            };
            results.push(new_item);
        }
        results
    }
}
