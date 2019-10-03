pub mod counts;
mod hashing;
pub mod minhashes;
pub mod scaled;

use needletail::SequenceRecord;

pub use hashing::ItemHash;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Hash, Serialize)]
pub struct KmerCount {
    pub hash: ItemHash,
    pub kmer: Vec<u8>,
    pub count: u64,
    pub extra_count: u64,
}

pub trait HashScheme {
    fn process(&mut self, seq: SequenceRecord);
    fn total_kmers(&self) -> usize;
    fn into_vec(self) -> Vec<KmerCount>;
}
