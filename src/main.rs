extern crate murmurhash3;
extern crate needletail;

mod minhash_queue;

use minhash_queue::MinHashQueue;

fn main() {
   let N_HASHES = 20000;
   let mut queue = MinHashQueue::new(N_HASHES);

}
