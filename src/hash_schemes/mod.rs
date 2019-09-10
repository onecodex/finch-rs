pub mod counts;
pub mod minhashes;
pub mod scaled;

use needletail::SequenceRecord;

// The individual items to store in the BinaryHeap
pub type ItemHash = u64;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Hash, Serialize)]
pub struct KmerCount {
    pub hash: ItemHash,
    pub kmer: Vec<u8>,
    pub count: u16,
    pub extra_count: u16,
}

pub trait HashScheme {
    fn process(&mut self, seq: SequenceRecord);
    fn total_kmers(&self) -> usize;
    fn into_vec(self) -> Vec<KmerCount>;
}
