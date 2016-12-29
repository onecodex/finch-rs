use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::usize;
use std::sync::Mutex;

use murmurhash3::murmurhash3_x64_128;
use needletail::bitkmer::{BitKmer, fast_extend_kmer};

fn to_ikmer(slice: &[u8]) -> BitKmer {
    let k = slice.len() as u8;

    let mut kmer = 0u64;
    for i in 0..k {
        kmer = fast_extend_kmer(&kmer, &k, &slice[i as usize]);
    }
    kmer
}

#[test]
fn test_to_ikmer() {
    let mut ikmer: BitKmer = to_ikmer("C".as_bytes());
    assert_eq!(1 as BitKmer, ikmer);

    ikmer = to_ikmer("TTA".as_bytes());
    assert_eq!(60 as BitKmer, ikmer);

    ikmer = to_ikmer("AAA".as_bytes());
    assert_eq!(0 as BitKmer, ikmer);
}


// TODO: add count to this?
#[derive(Debug)]
pub struct HashedItem<T> {
    hash: usize,
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
        other.hash.cmp(&self.hash)
    }
}

impl<T> PartialOrd for HashedItem<T> {
    fn partial_cmp(&self, other: &HashedItem<T>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


pub struct MinHashQueue {
    heap: BinaryHeap<HashedItem<BitKmer>>,
    size: usize,
    lock: Mutex<()>,
}

impl MinHashQueue {
    pub fn new(size: usize) -> Self {
        MinHashQueue {
            heap: BinaryHeap::new(),
            size: size,
            lock: Mutex::new(()),
        }
    }

    pub fn push(&mut self, data: &[u8]) {
        let new_hash = HashedItem {
            hash: murmurhash3_x64_128(data, 42).0 as usize,
            item: to_ikmer(data),
        };
        let add_hash = match self.heap.peek() {
            None => true,
            Some(old_max_hash) => new_hash < *old_max_hash,
        };
        if add_hash {
            self.lock.lock().unwrap();
            self.heap.push(new_hash);
            if self.heap.len() > self.size {
                self.heap.pop();
            }
        }
    }

    pub fn into_vec(self) -> Vec<HashedItem<BitKmer>> {
        let mut vec = self.heap.into_sorted_vec();
        vec.reverse();
        vec
    }
}

#[test]
fn test_queue() {
    let mut queue = MinHashQueue::new(2);
    queue.push(b"ca");
    queue.push(b"cc");
    queue.push(b"ac");
    let array = queue.into_vec();
    assert_eq!(array[0].item, 4u64);
    assert_eq!(array[1].item, 1u64);
}
