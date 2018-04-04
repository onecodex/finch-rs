use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::hash::{Hasher, BuildHasherDefault};
use std::usize;

use murmurhash3::murmurhash3_x64_128;


// The individual items to store in the BinaryHeap
pub type ItemHash = usize;

#[inline]
pub fn hash_f(item: &[u8], seed: u64) -> ItemHash {
    murmurhash3_x64_128(item, seed).0 as ItemHash
}


#[derive(Debug, Clone)]
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


#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Hash, Serialize)]
pub struct KmerCount {
    pub hash: ItemHash,
    pub kmer: Vec<u8>,
    pub count: u16,
    pub extra_count: u16,
}


#[derive(Clone)]
pub struct MinHashKmers {
    hashes: BinaryHeap<HashedItem<Vec<u8>>>,
    counts: HashMap<ItemHash, (u16, u16), BuildHasherDefault<NoHashHasher>>,
    total_count: u64,
    size: usize,
    seed: u64,
    // heap_lock: Mutex<()>,
    // instead of map_lock, look into using https://docs.rs/chashmap/0.1.2/chashmap/
    // map_lock: Mutex<()>,
}

impl MinHashKmers {
    pub fn new(size: usize, seed: u64) -> Self {
        MinHashKmers {
            hashes: BinaryHeap::with_capacity(size + 1),
            counts: HashMap::with_capacity_and_hasher(size, BuildHasherDefault::default()),
            total_count: 0,
            size: size,
            seed: seed,
            // heap_lock: Mutex::new(()),
            // map_lock: Mutex::new(()),
        }
    }

    pub fn push(&mut self, kmer: &[u8], extra_count: u8) {
        let new_hash = hash_f(kmer, self.seed);
        let add_hash = match self.hashes.peek() {
            None => true,
            Some(old_max_hash) => (new_hash <= (*old_max_hash).hash) || (self.hashes.len() < self.size),
        };

        if add_hash {
            self.total_count += 1;
            if self.counts.contains_key(&new_hash) {
                // let _lock = self.map_lock.lock().unwrap();
                let count = self.counts.entry(new_hash).or_insert((0u16, 0u16));
                (*count).0 += 1;
                (*count).1 += extra_count as u16;
                // drop(_lock);
            } else {
                // let _ = self.heap_lock.lock().unwrap();
                self.hashes.push(HashedItem {
                    hash: new_hash,
                    item: kmer.to_owned(),
                });
                // let _map_lock = self.map_lock.lock().unwrap();
                self.counts.insert(new_hash, (1u16, extra_count as u16));
                if self.hashes.len() > self.size {
                    let hash = self.hashes.pop().unwrap();
                    let old_count = self.counts.remove(&hash.hash).unwrap().0;
                    // TODO: check if old_count / total_count > ...?
                    self.total_count -= old_count as u64;
                }
                // drop(_lock);
                // drop(_map_lock);
            }
        }
    }

    pub fn total_kmers(&self) -> usize {
        self.total_count as usize
    }

    pub fn into_vec(self) -> Vec<KmerCount> {
        let mut vec = self.hashes.into_sorted_vec();

        let mut results = Vec::with_capacity(vec.len());
        for item in vec.drain(..) {
            let counts = *self.counts.get(&item.hash).unwrap();
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

#[test]
fn test_minhashkmers() {
    let mut queue = MinHashKmers::new(3, 42);
    queue.push(b"ca", 0);
    queue.push(b"cc", 1);
    queue.push(b"ac", 0);
    queue.push(b"ac", 1);
    let array = queue.into_vec();
    assert_eq!(array[0].kmer, vec![b'c', b'c']);
    assert_eq!(array[0].count, 1u16);
    assert_eq!(array[0].extra_count, 1u16);
    assert!(array[0].hash < array[1].hash);
    assert_eq!(array[1].kmer, vec![b'c', b'a']);
    assert_eq!(array[1].count, 1u16);
    assert_eq!(array[1].extra_count, 0u16);
    assert!(array[1].hash < array[2].hash);
    assert_eq!(array[2].kmer, vec![b'a', b'c']);
    assert_eq!(array[2].count, 2u16);
    assert_eq!(array[2].extra_count, 1u16);
}

#[test]
fn test_longer_sequence() {
    let mut queue = MinHashKmers::new(100, 42);

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
