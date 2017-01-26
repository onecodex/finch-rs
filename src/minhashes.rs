use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::hash::{Hasher, BuildHasherDefault};
use std::usize;
use std::sync::Mutex;

use murmurhash3::murmurhash3_x64_128;
use needletail::bitkmer::{BitKmer, str_to_bitmer};


// The individual items to store in the BinaryHeap
pub type ItemHash = usize;

#[inline]
pub fn hash_f(item: &[u8]) -> ItemHash {
    murmurhash3_x64_128(item, 42).0 as ItemHash
}


#[derive(Debug)]
struct HashedItem<T> {
    hash: ItemHash,
    item: T,
}


impl<T> PartialEq for HashedItem<T> {
    fn eq(&self, other: &HashedItem<T>) -> bool {
        other.hash.eq(&self.hash)
    }
}

impl<T> Eq for HashedItem<T> {}

impl<T> Ord for HashedItem<T> {
    fn cmp(&self, other: &HashedItem<T>) -> Ordering {
        self.hash.cmp(&other.hash)
    }
}

impl<T> PartialOrd for HashedItem<T> {
    fn partial_cmp(&self, other: &HashedItem<T>) -> Option<Ordering> {
        Some(self.hash.cmp(&other.hash))
    }
}


/// If we're using a `HashMap` where the keys themselves are hashes, it's
/// a little silly to re-hash them. That's where the `NoHashHasher` comes in.
pub struct NoHashHasher(u64);

impl Default for NoHashHasher {
    #[inline]
    fn default() -> NoHashHasher {
        NoHashHasher(0x0000000000000000)
    }
}

impl Hasher for NoHashHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        *self = NoHashHasher(
            ((bytes[0] as u64) << 24) +
            ((bytes[1] as u64) << 16) +
            ((bytes[2] as u64) << 8) +
            (bytes[3] as u64)
        );
    }
    fn finish(&self) -> u64 { self.0 }
}


#[derive(Debug)]
pub struct KmerCount {
    pub hash: ItemHash,
    pub kmer: BitKmer,
    pub count: u16,
}


pub struct MinHashKmers {
    hashes: BinaryHeap<HashedItem<BitKmer>>,
    counts: HashMap<ItemHash, u16, BuildHasherDefault<NoHashHasher>>,
    size: usize,
    heap_lock: Mutex<()>,
    map_lock: Mutex<()>,
}

impl MinHashKmers {
    pub fn new(size: usize) -> Self {
        MinHashKmers {
            hashes: BinaryHeap::with_capacity(size + 1),
            counts: HashMap::with_capacity_and_hasher(size, BuildHasherDefault::default()),
            size: size,
            heap_lock: Mutex::new(()),
            map_lock: Mutex::new(()),
        }
    }

    pub fn push(&mut self, kmer: &[u8]) {
        let new_hash = hash_f(kmer);
        let add_hash = match self.hashes.peek() {
            None => true,
            Some(old_max_hash) => (new_hash <= (*old_max_hash).hash) || (self.hashes.len() < self.size),
        };

        if add_hash {
            let new_hash_item = HashedItem {
                hash: new_hash,
                item: str_to_bitmer(kmer),
            };
            if self.counts.contains_key(&new_hash) {
                self.map_lock.lock().unwrap();
                let count = self.counts.entry(new_hash).or_insert(0u16);
                *count += 1;
            } else {
                self.heap_lock.lock().unwrap();
                self.hashes.push(new_hash_item);
                self.counts.insert(new_hash, 1u16);
                if self.hashes.len() > self.size {
                    let hash = self.hashes.pop().unwrap();
                    self.counts.remove(&hash.hash);
                }
            }
        }
    }

    pub fn into_vec(self) -> Vec<KmerCount> {
        let mut vec = self.hashes.into_sorted_vec();

        let mut results = Vec::with_capacity(vec.len());
        for item in &vec {
            let new_item = KmerCount {
                hash: item.hash,
                kmer: item.item,
                count: *self.counts.get(&item.hash).unwrap(),
            };
            results.push(new_item);
        }
        results
    }
}

#[test]
fn test_minhashkmers() {
    let mut queue = MinHashKmers::new(2);
    queue.push(b"ca");
    queue.push(b"cc");
    queue.push(b"ac");
    queue.push(b"ac");
    let array = queue.into_vec();
    assert_eq!(array[0].kmer.0, 4u64);
    assert_eq!(array[0].count, 1u16);
    assert_eq!(array[1].kmer.0, 1u64);
    assert_eq!(array[1].count, 2u16);
}

#[test]
fn test_longer_sequence() {
    let mut queue = MinHashKmers::new(100);

    // for "ACACGGAAATCCTCACGTCGCGGCGCCGGGC"

    // hashes should be:
    //     (3186265289206375993,
    //      3197567229193635484,
    //      5157287830980272133,
    //      7515070071080094037,
    //      9123665698461883699,
    //      9650810550987401968,
    //      10462414310441547028,
    //      12872951831549606632,
    //      13584836512372089324,
    //      14093285637546356047,
    //      16069721578136260683)

}
