use std::cmp::Ordering;
use std::hash::Hasher;

use murmurhash3::murmurhash3_x64_128;

// The individual items to store in the BinaryHeap
pub type ItemHash = u64;

#[inline]
pub fn hash_f(item: &[u8], seed: u64) -> ItemHash {
    murmurhash3_x64_128(item, seed).0
}

#[derive(Debug, Clone)]
pub(crate) struct HashedItem<T> {
    pub(crate) hash: ItemHash,
    pub(crate) item: T,
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
#[derive(Default)]
pub struct NoHashHasher(u64);

impl Hasher for NoHashHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        *self = NoHashHasher(
            (u64::from(bytes[0]) << 24)
                + (u64::from(bytes[1]) << 16)
                + (u64::from(bytes[2]) << 8)
                + u64::from(bytes[3]),
        );
    }
    fn finish(&self) -> u64 {
        self.0
    }
}
