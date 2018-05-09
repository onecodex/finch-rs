use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::io::Result as IOResult;

#[cfg(feature = "mash_format")]
use capnp::{serialize_packed, message};
use needletail::bitkmer::{bytes_to_bitmer, bitmer_to_bytes};
use serde::de::{Deserialize, Deserializer, Visitor, SeqAccess, MapAccess};
use serde::ser::{Serialize, Serializer, SerializeStruct};

use filtering::{FilterParams, filter_sketch};
use minhashes::{KmerCount, hash_f};
#[cfg(feature = "mash_format")]
use mash_capnp::{min_hash};
use ::Result as FinchResult;
// min_hash, reference_list, reference, locus_list, locus


#[allow(non_snake_case)]
#[derive(Debug, Serialize, Deserialize)]
pub struct SketchDistance {
    pub containment: f64,
    pub jaccard: f64,
    pub mashDistance: f64,
    pub commonHashes: u64,
    pub totalHashes: u64,
    pub query: String,
    pub reference: String,
}


#[allow(non_snake_case)]
#[derive(Debug, Eq, PartialEq, Clone )]
pub struct Sketch {
    pub name: String,
    pub seqLength: Option<u64>,
    pub numValidKmers: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    hashes: Vec<KmerCount>,
}


impl Serialize for Sketch {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut hash_list = Vec::with_capacity(self.hashes.len());
        let mut kmer_list = Vec::with_capacity(self.hashes.len());
        let mut count_list = Vec::with_capacity(self.hashes.len());
        for hash in &self.hashes {
            hash_list.push(hash.hash.to_string());
            kmer_list.push(String::from_utf8(hash.kmer.clone()).unwrap());
            count_list.push(hash.count);
        }

        let mut state = serializer.serialize_struct("Sketch", 8)?;
        state.serialize_field("name", &self.name)?;
        state.serialize_field("seqLength", &self.seqLength)?;
        state.serialize_field("numValidKmers", &self.numValidKmers)?;
        state.serialize_field("comment", &self.comment)?;
        state.serialize_field("filters", &self.filters)?;
        state.serialize_field("hashes", &hash_list)?;
        state.serialize_field("kmers", &kmer_list)?;
        state.serialize_field("counts", &count_list)?;
        state.end()
    }
}


#[allow(non_snake_case)]
#[derive(Debug, Deserialize)]
struct JSONSketch {
    pub name: String,
    pub seqLength: Option<u64>,
    pub numValidKmers: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    hashes: Vec<String>,
    kmers: Option<Vec<String>>,
    counts: Option<Vec<u16>>,
}


impl<'de> Deserialize<'de> for Sketch {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>
    {
        let jsketch = JSONSketch::deserialize(deserializer)?;

        let mut kmercount_list = Vec::with_capacity(jsketch.hashes.len());
        for i in 0..jsketch.hashes.len() {
            let hash;
            match jsketch.hashes[i].parse::<usize>() {
                Ok(t) => hash = t,
                Err(_) => break,
            }
            let kmer;
            match jsketch.kmers {
                Some(ref v) => kmer = v[i].clone().into_bytes(),
                None => kmer = Vec::new(), // TODO: this might be slow?
            }
            let count;
            match jsketch.counts {
                Some(ref v) => count = v[i],
                None => count = 1,
            }
            kmercount_list.push(KmerCount {
                hash: hash,
                kmer: kmer,
                count: count,
                extra_count: count / 2,
            });
        }
        Ok(Sketch {
            name: jsketch.name,
            seqLength: jsketch.seqLength,
            numValidKmers: jsketch.numValidKmers,
            comment: jsketch.comment,
            filters: jsketch.filters,
            hashes: kmercount_list,
        })
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Serialize, Deserialize)]
pub struct MultiSketch {
    pub kmer: u8,
    pub alphabet: String,
    pub preserveCase: bool,
    pub canonical: bool,
    pub sketchSize: u32,
    pub hashType: String,
    pub hashBits: u16,
    pub hashSeed: u64,
    pub sketches: Vec<Sketch>,
}

impl Sketch {
    pub fn new(name: &str, length: u64, n_kmers: u64, kmercounts: Vec<KmerCount>, filters: &HashMap<String, String>) -> Self {
        Sketch {
            name: String::from(name),
            seqLength: Some(length),
            numValidKmers: Some(n_kmers),
            comment: Some(String::from("")),
            filters: Some(filters.clone()),
            hashes: kmercounts,
        }
    }

    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn get_kmers(&self) -> Option<Vec<KmerCount>> {
        Some(self.hashes.clone())
    }

    pub fn apply_filtering(&mut self, filters: &FilterParams) -> bool {
        let (filtered_hashes, filter_stats) = filter_sketch(&self.hashes, &filters);
        self.hashes = filtered_hashes;
        self.filters = Some(filter_stats);
        true
    }
}
