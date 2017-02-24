use std::collections::HashMap;
use needletail::bitkmer::{str_to_bitmer, bitmer_to_str};

use minhashes::{KmerCount, hash_f};

#[allow(non_snake_case)]
#[derive(Serialize, Deserialize)]
pub struct SketchDistance {
    pub containment: f64,
    pub jaccard: f64,
    pub mashDistance: f64,
    pub common_hashes: u64,
    pub total_hashes: u64,
    pub query: String,
    pub reference: String,
}

#[allow(non_snake_case)]
#[derive(Serialize, Deserialize)]
pub struct JSONMultiSketch {
    pub kmer: u8,
    pub alphabet: String,
    pub preserveCase: bool,
    pub canonical: bool,
    pub sketchSize: u32,
    pub hashType: String,
    pub hashBits: u16,
    pub hashSeed: u64,
    pub sketches: Vec<JSONSketch>,
}

#[allow(non_snake_case)]
#[derive(Deserialize, Eq, PartialEq, Clone, Serialize)]
pub struct JSONSketch {
    pub name: String,
    seqLength: Option<u64>,
    comment: Option<String>,
    filters: Option<HashMap<String, String>>,
    hashes: Vec<String>,
    kmers: Option<Vec<String>>,
    counts: Option<Vec<u16>>,
}

impl JSONSketch {
    pub fn new(name: &str, length: u64, kmercounts: Vec<KmerCount>, filters: &HashMap<String, String>) -> Self {
        let mut hash_list = Vec::with_capacity(kmercounts.len());
        let mut kmer_list = Vec::with_capacity(kmercounts.len());
        let mut count_list = Vec::with_capacity(kmercounts.len());
        for hash in &kmercounts {
            hash_list.push(hash.hash.to_string());
            kmer_list.push(String::from_utf8(hash.kmer.clone()).unwrap());
            count_list.push(hash.count);
        }
        JSONSketch {
            name: String::from(name),
            seqLength: Some(length),
            comment: Some(String::from("")),
            filters: Some(filters.clone()),
            hashes: hash_list,
            kmers: Some(kmer_list),
            counts: Some(count_list),
        }
    }

    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn get_kmers(&self) -> Option<Vec<KmerCount>> {
        let mut kmercount_list = Vec::with_capacity(self.hashes.len());
        for i in 0..self.hashes.len() {
            let hash;
            match self.hashes[i].parse::<usize>() {
                Ok(t) => hash = t,
                Err(_) => return None,
            }
            let kmer;
            match self.kmers {
                Some(ref v) => kmer = v[i].clone().into_bytes(),
                None => return None,
            }
            let count;
            match self.counts {
                Some(ref v) => count = v[i],
                None => return None,
            }
            kmercount_list.push(KmerCount {
                hash: hash,
                kmer: kmer,
                count: count,
            });
        }
        Some(kmercount_list)
    }
}


#[repr(C, packed)]
pub struct BinarySketch {
    len: u32,
    kmer_size: u8,
    kmers: Box<[u64]>,
    counts: Box<[u16]>,
}

impl BinarySketch {
    pub fn new(kmercounts: Vec<KmerCount>) -> Self {
        let mut kmer_list = Vec::with_capacity(kmercounts.len());
        let mut count_list = Vec::with_capacity(kmercounts.len());
        for hash in &kmercounts {
            kmer_list.push(str_to_bitmer(&hash.kmer).0 as u64);
            count_list.push(hash.count);
        }
        BinarySketch {
            len: kmercounts.len() as u32,
            kmer_size: kmercounts[0].kmer.len() as u8,
            kmers: kmer_list.into_boxed_slice(), 
            counts: count_list.into_boxed_slice(),
        }
    }

    pub fn get_kmers(&self) -> Option<Vec<KmerCount>> {
        let mut kmercounts = Vec::with_capacity(self.len as usize);
        for i in 0..self.len {
            let bitmer = (*self.kmers)[i as usize];
            let kmer = bitmer_to_str((bitmer, self.kmer_size));
            kmercounts.push(KmerCount {
                // there's an assumption here that the seed is 42
                hash: hash_f(&kmer, 42),
                kmer: kmer,
                count: (*self.counts)[i as usize],
            });
        }
        Some(kmercounts)
    }
}
