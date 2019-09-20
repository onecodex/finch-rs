use std::collections::{BinaryHeap, HashMap};
use std::hash::BuildHasherDefault;
use std::usize;

use needletail::{Sequence, SequenceRecord};

use crate::hash_schemes::minhashes::{hash_f, HashedItem, MinHashKmers, NoHashHasher};
use crate::hash_schemes::{HashScheme, ItemHash, KmerCount};

#[derive(Clone, Debug)]
pub struct ScaledKmers {
    hashes: BinaryHeap<HashedItem<Vec<u8>>>,
    counts: HashMap<ItemHash, (u16, u16), BuildHasherDefault<NoHashHasher>>,
    kmer_length: u8,
    total_kmers: u64,
    size: usize,
    max_hash: u64,
    seed: u64,
}

impl ScaledKmers {
    pub fn new(size: usize, scale: f64, kmer_length: u8, seed: u64) -> Self {
        let iscale = (1. / scale) as u64;
        ScaledKmers {
            hashes: BinaryHeap::with_capacity(size),
            counts: HashMap::with_capacity_and_hasher(size, BuildHasherDefault::default()),
            kmer_length,
            total_kmers: 0,
            size,
            max_hash: u64::max_value() / iscale,
            seed,
        }
    }

    #[allow(clippy::map_entry)]
    pub fn push(&mut self, kmer: &[u8], extra_count: u8) {
        self.total_kmers += 1;
        let new_hash = hash_f(kmer, self.seed);

        if new_hash <= self.max_hash || (self.hashes.len() <= self.size && self.size != 0) {
            if self.counts.contains_key(&new_hash) {
                let count = self.counts.entry(new_hash).or_insert((0u16, 0u16));
                (*count).0 += 1;
                (*count).1 += u16::from(extra_count);
            } else {
                self.hashes.push(HashedItem {
                    hash: new_hash,
                    item: kmer.to_owned(),
                });
                self.counts.insert(new_hash, (1u16, u16::from(extra_count)));
                if self.hashes.len() > self.size
                    && (*self.hashes.peek().unwrap()).hash > self.max_hash
                {
                    let hash = self.hashes.pop().unwrap();
                    let _ = self.counts.remove(&hash.hash).unwrap();
                }
            }
        }
    }
}

impl HashScheme for ScaledKmers {
    fn process(&mut self, seq: SequenceRecord) {
        let rc = seq.reverse_complement();
        for (_, kmer, is_rev_complement) in
            seq.normalize(false).canonical_kmers(self.kmer_length, &rc)
        {
            let rc_count = if is_rev_complement { 1u8 } else { 0u8 };
            self.push(kmer, rc_count);
        }
    }

    fn total_kmers(&self) -> usize {
        self.total_kmers as usize
    }

    fn into_vec(self) -> Vec<KmerCount> {
        let mut vec = self.hashes.into_sorted_vec();

        let mut results = Vec::with_capacity(vec.len());
        for item in vec.drain(..) {
            let counts = self.counts[&item.hash];
            let new_item = KmerCount {
                hash: item.hash,
                kmer: item.item,
                count: counts.0,
                extra_count: counts.1,
            };
            results.push(new_item);
        }
        results
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
            };
            results.push(new_item);
        }
        results
    }
}

impl From<ScaledKmers> for MinHashKmers {
    fn from(value: ScaledKmers) -> Self {
        let size = value.size;
        let kmer_length = value.kmer_length;
        let total_kmers = value.total_kmers;
        let seed = value.seed;

        let mut hashes = value.hashes.into_sorted_vec();
        hashes.truncate(size);

        let mut counts = value.counts;
        counts.retain(|&k, _| hashes.binary_search_by_key(&k, |a| a.hash).is_ok());

        MinHashKmers {
            hashes: hashes.into(),
            counts,
            kmer_length,
            total_kmers,
            size,
            seed,
        }
    }
}

#[cfg(test)]
mod test {
    use proptest::prelude::*;

    use crate::distance::distance_scaled;

    use super::*;

    #[test]
    fn test_minhashkmers_scaled_1() {
        // Scaled=1 should hold all possible kmers
        let mut queue = ScaledKmers::new(3, 1., 2, 42);
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.into_vec();
        assert_eq!(array[0].kmer, b"cc");
        assert_eq!(array[0].count, 1u16);
        assert_eq!(array[0].extra_count, 1u16);
        assert!(array[0].hash < array[1].hash);
        assert_eq!(array[1].kmer, b"ca");
        assert_eq!(array[1].count, 1u16);
        assert_eq!(array[1].extra_count, 0u16);
        assert!(array[1].hash < array[2].hash);
        assert_eq!(array[2].kmer, b"ac");
        assert_eq!(array[2].count, 2u16);
        assert_eq!(array[2].extra_count, 1u16);
    }

    #[test]
    fn test_minhashkmers_scaled_1000() {
        // Scaled=1000 should exclude all these hashes,
        // but since only 3 are added and size==3 they should all be present
        let mut queue = ScaledKmers::new(3, 0.001, 2, 42);
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.into_vec();
        assert_eq!(array[0].kmer, b"cc");
        assert_eq!(array[0].count, 1u16);
        assert_eq!(array[0].extra_count, 1u16);
        assert!(array[0].hash < array[1].hash);
        assert_eq!(array[1].kmer, b"ca");
        assert_eq!(array[1].count, 1u16);
        assert_eq!(array[1].extra_count, 0u16);
        assert!(array[1].hash < array[2].hash);
        assert_eq!(array[2].kmer, b"ac");
        assert_eq!(array[2].count, 2u16);
        assert_eq!(array[2].extra_count, 1u16);
    }

    #[test]
    fn test_minhashkmers_eviction() {
        // try again, but evict one of the kmers
        let mut queue = ScaledKmers::new(1, 0.01, 4, 42);
        // random kmer that hashes above max_hash
        queue.push(b"AAAA", 0);
        // now fill with kmers that hash below to evict it
        queue.push(b"AGTA", 0);
        queue.push(b"CCCC", 1);
        queue.push(b"ATAA", 0);
        let array = queue.into_vec();
        assert_eq!(array.len(), 3, "Only small hashes should be left");
        assert!(array.iter().all(|e| e.kmer != b"AAAA"))
    }

    #[test]
    fn test_minhashkmers_eviction_and_conversion() {
        let mut queue = ScaledKmers::new(4, 0.01, 4, 42);
        // random kmer that hashes above max_hash
        queue.push(b"AAAA", 0);
        // now fill with kmers that hash below max_hash.
        queue.push(b"AGTA", 0);
        queue.push(b"CCCC", 1);
        queue.push(b"ATAA", 0);

        // The scaled minhash should have 3 hashes, but since we asked for
        // size=4 MinHashKmers should have size=4
        let mh: MinHashKmers = queue.into();
        let array = mh.into_vec();

        assert_eq!(array.len(), 4, "Should keep all four hashes");
    }

    #[test]
    fn test_minhashkmers_pure_scaled_empty() {
        let mut queue = ScaledKmers::new(0, 0.001, 2, 42);
        // all these hashes are out of range for scaled=1000
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.into_vec();
        assert_eq!(array.len(), 0);
    }

    #[test]
    fn test_minhashkmers_pure_scaled() {
        let mut queue = ScaledKmers::new(0, 0.001, 2, 42);
        // all these hashes are out of range
        queue.push(b"ca", 0);
        queue.push(b"cc", 1);
        queue.push(b"ac", 0);
        queue.push(b"ac", 1);
        let array = queue.into_vec();
        assert_eq!(array.len(), 0);
    }

    proptest! {
        #[test]
        fn pure_scaled_check(seq in "[ACGT]{500,}") {
            let mut queue = ScaledKmers::new(0, 1. / 100., 2, 42);
            let max_hash = u64::max_value() / 100;
            for kmer in seq.as_bytes().windows(4) {
                queue.push(kmer, 0);
            }
            let array = queue.into_vec();
            assert!(array.iter().all(|item| item.hash <= max_hash));
        }
    }

    #[test]
    fn test_distance_scaled() -> Result<(), Box<dyn std::error::Error>> {
        let mut queue1 = ScaledKmers::new(3, 0.001, 2, 42);
        queue1.push(b"ca", 0);
        queue1.push(b"cc", 1);
        queue1.push(b"ac", 0);
        queue1.push(b"ac", 1);
        let array1 = queue1.into_vec();

        let mut queue2 = ScaledKmers::new(3, 0.001, 2, 42);
        queue2.push(b"ca", 0);
        queue2.push(b"cc", 1);
        queue2.push(b"ac", 0);
        queue2.push(b"ac", 1);
        let array2 = queue2.into_vec();

        let dist = distance_scaled(&array1, &array2, "", "")?;
        assert_eq!(dist.jaccard, 1.0);
        assert_eq!(dist.containment, 1.0);
        assert_eq!(dist.commonHashes, 3);

        Ok(())
    }
}
