use std::collections::{BinaryHeap, HashMap};
use std::hash::BuildHasherDefault;

use needletail::Sequence;

use crate::sketch_schemes::hashing::{hash_f, HashedItem, NoHashHasher};
use crate::sketch_schemes::{ItemHash, KmerCount, SketchParams, SketchScheme};

#[derive(Clone, Debug)]
pub struct MashSketcher {
    hashes: BinaryHeap<HashedItem<Vec<u8>>>,
    counts: HashMap<ItemHash, (u32, u32), BuildHasherDefault<NoHashHasher>>,
    kmer_length: u8,
    total_kmers: u64,
    total_bases: u64,
    size: usize,
    seed: u64,
}

impl MashSketcher {
    pub fn new(size: usize, kmer_length: u8, seed: u64) -> Self {
        MashSketcher {
            hashes: BinaryHeap::with_capacity(size + 1),
            counts: HashMap::with_capacity_and_hasher(size, BuildHasherDefault::default()),
            kmer_length,
            total_kmers: 0,
            total_bases: 0,
            size,
            seed,
        }
    }

    #[allow(clippy::map_entry)]
    pub fn push(&mut self, kmer: &[u8], extra_count: u8) {
        self.total_kmers += 1;
        let new_hash = hash_f(kmer, self.seed);
        let add_hash = match self.hashes.peek() {
            None => true,
            Some(old_max_hash) => {
                (new_hash <= old_max_hash.hash) || (self.hashes.len() < self.size)
            }
        };

        if add_hash {
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
                if self.hashes.len() > self.size {
                    let hash = self.hashes.pop().unwrap();
                    let _ = self.counts.remove(&hash.hash).unwrap();
                }
            }
        }
    }
}

impl SketchScheme for MashSketcher {
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
        SketchParams::Mash {
            kmers_to_sketch: self.size,
            final_size: self.size,
            no_strict: false,
            kmer_length: self.kmer_length,
            hash_seed: self.seed,
        }
    }
}

#[test]
fn test_minhashkmers() {
    let mut queue = MashSketcher::new(3, 2, 42);
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

//#[test]
//fn test_longer_sequence() {
//    let mut queue = MinHashKmers::new(100, 21, 42);
//
//    // for "ACACGGAAATCCTCACGTCGCGGCGCCGGGC"
//
//    // hashes should be:
//    //     (3186265289206375993,
//    //      3197567229193635484,
//    //      5157287830980272133,
//    //      7515070071080094037,
//    //      9123665698461883699,
//    //      9650810550987401968,
//    //      10462414310441547028,
//    //      12872951831549606632,
//    //      13584836512372089324,
//    //      14093285637546356047,
//    //      16069721578136260683)
//}
