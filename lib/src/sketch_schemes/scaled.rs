use std::collections::{BinaryHeap, HashMap};
use std::hash::BuildHasherDefault;

use needletail::Sequence;

use crate::sketch_schemes::hashing::{hash_f, HashedItem, NoHashHasher};
use crate::sketch_schemes::{ItemHash, KmerCount, SketchParams, SketchScheme};

#[derive(Clone, Debug)]
pub struct ScaledSketcher {
    hashes: BinaryHeap<HashedItem<Vec<u8>>>,
    counts: HashMap<ItemHash, (u32, u32), BuildHasherDefault<NoHashHasher>>,
    kmer_length: u8,
    total_kmers: u64,
    total_bases: u64,
    size: usize,
    max_hash: u64,
    seed: u64,
}

impl ScaledSketcher {
    pub fn new(size: usize, scale: f64, kmer_length: u8, seed: u64) -> Self {
        let iscale = (1. / scale) as u64;
        ScaledSketcher {
            hashes: BinaryHeap::with_capacity(size),
            counts: HashMap::with_capacity_and_hasher(size, BuildHasherDefault::default()),
            kmer_length,
            total_kmers: 0,
            total_bases: 0,
            size,
            max_hash: u64::MAX / iscale,
            seed,
        }
    }

    #[allow(clippy::map_entry)]
    pub fn push(&mut self, kmer: &[u8], extra_count: u8) {
        self.total_kmers += 1;
        let new_hash = hash_f(kmer, self.seed);

        if new_hash <= self.max_hash || (self.hashes.len() <= self.size && self.size != 0) {
            if self.counts.contains_key(&new_hash) {
                let count = self.counts.entry(new_hash).or_insert((0, 0));
                (*count) = (
                    count.0.saturating_add(1),
                    count.1.saturating_add(u32::from(extra_count)),
                );
            } else {
                self.hashes.push(HashedItem {
                    hash: new_hash,
                    item: kmer.to_owned(),
                });
                self.counts.insert(new_hash, (1, u32::from(extra_count)));
                if self.hashes.len() > self.size && self.hashes.peek().unwrap().hash > self.max_hash
                {
                    let hash = self.hashes.pop().unwrap();
                    let _ = self.counts.remove(&hash.hash).unwrap();
                }
            }
        }
    }
}

impl SketchScheme for ScaledSketcher {
    fn process<'seq, 'a, 'inner>(&'a mut self, seq: &'seq dyn Sequence<'inner>)
    where
        'a: 'seq,
        'seq: 'inner,
    {
        self.total_bases += seq.sequence().len() as u64;
        let norm_seq = seq.normalize(false);

        let rc = norm_seq.reverse_complement();
        for (_, kmer, is_rev_complement) in norm_seq.canonical_kmers(self.kmer_length, &rc) {
            let rc_count = u8::from(is_rev_complement);
            self.push(kmer, rc_count);
        }
    }

    fn total_bases_and_kmers(&self) -> (u64, u64) {
        (self.total_bases, self.total_kmers)
    }

    fn to_vec(&self) -> Vec<KmerCount> {
        let mut vec = self.hashes.clone().into_sorted_vec();

        let mut results = Vec::with_capacity(vec.len());
        for item in vec.drain(..) {
            let counts = self.counts[&item.hash];
            let new_item = KmerCount {
                hash: item.hash,
                kmer: item.item,
                count: counts.0,
                extra_count: counts.1,
                label: None,
            };
            results.push(new_item);
        }
        results
    }

    fn parameters(&self) -> SketchParams {
        SketchParams::Scaled {
            kmers_to_sketch: self.size,
            kmer_length: self.kmer_length,
            scale: 1. / (u64::MAX as f64 / self.max_hash as f64),
            hash_seed: self.seed,
        }
    }
}

#[cfg(test)]
mod test {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_minhashkmers_scaled_1() {
        // Scaled=1 should hold all possible kmers
        let mut queue = ScaledSketcher::new(3, 1., 2, 42);
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.to_vec();
        assert_eq!(array[0].kmer, b"cc");
        assert_eq!(array[0].count, 1);
        assert_eq!(array[0].extra_count, 1);
        assert!(array[0].hash < array[1].hash);
        assert_eq!(array[1].kmer, b"ca");
        assert_eq!(array[1].count, 1);
        assert_eq!(array[1].extra_count, 0);
        assert!(array[1].hash < array[2].hash);
        assert_eq!(array[2].kmer, b"ac");
        assert_eq!(array[2].count, 2);
        assert_eq!(array[2].extra_count, 1);
    }

    #[test]
    fn test_minhashkmers_scaled_1000() {
        // Scaled=1000 should exclude all these hashes,
        // but since only 3 are added and size==3 they should all be present
        let mut queue = ScaledSketcher::new(3, 0.001, 2, 42);
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.to_vec();
        assert_eq!(array[0].kmer, b"cc");
        assert_eq!(array[0].count, 1);
        assert_eq!(array[0].extra_count, 1);
        assert!(array[0].hash < array[1].hash);
        assert_eq!(array[1].kmer, b"ca");
        assert_eq!(array[1].count, 1);
        assert_eq!(array[1].extra_count, 0);
        assert!(array[1].hash < array[2].hash);
        assert_eq!(array[2].kmer, b"ac");
        assert_eq!(array[2].count, 2);
        assert_eq!(array[2].extra_count, 1);
    }

    #[test]
    fn test_minhashkmers_eviction() {
        // try again, but evict one of the kmers
        let mut queue = ScaledSketcher::new(1, 0.01, 4, 42);
        // random kmer that hashes above max_hash
        queue.push(b"AAAA", 0);
        // now fill with kmers that hash below to evict it
        queue.push(b"AGTA", 0);
        queue.push(b"CCCC", 1);
        queue.push(b"ATAA", 0);
        let array = queue.to_vec();
        assert_eq!(array.len(), 3, "Only small hashes should be left");
        assert!(array.iter().all(|e| e.kmer != b"AAAA"))
    }

    #[test]
    fn test_minhashkmers_pure_scaled_empty() {
        let mut queue = ScaledSketcher::new(0, 0.001, 2, 42);
        // all these hashes are out of range for scaled=1000
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.to_vec();
        assert_eq!(array.len(), 0);
    }

    #[test]
    fn test_minhashkmers_pure_scaled() {
        let mut queue = ScaledSketcher::new(0, 0.001, 2, 42);
        // all these hashes are out of range
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.to_vec();
        assert_eq!(array.len(), 0);
    }

    proptest! {
        #[test]
        fn pure_scaled_check(seq in "[ACGT]{500,}") {
            let mut queue = ScaledSketcher::new(0, 1. / 100., 2, 42);
            let max_hash = u64::max_value() / 100;
            for kmer in seq.as_bytes().windows(4) {
                queue.push(kmer, 0);
            }
            let array = queue.to_vec();
            assert!(array.iter().all(|item| item.hash <= max_hash));
        }
    }
}
