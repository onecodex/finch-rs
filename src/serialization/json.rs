use std::collections::HashMap;
use std::fmt;
use std::mem;

use failure::bail;
use serde::de::{self, Deserialize, Deserializer, Visitor};
use serde::ser::{Serialize, SerializeStruct, Serializer};

use crate::filtering::FilterParams;
pub use crate::serialization::mash::{read_mash_file, write_mash_file};
use crate::sketch_schemes::{KmerCount, SketchParams};
use crate::Result as FinchResult;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sketch {
    pub name: String,
    pub seq_length: Option<u64>,
    pub num_valid_kmers: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    pub hashes: Vec<KmerCount>,
}

impl Sketch {
    pub fn new(
        name: &str,
        length: u64,
        n_kmers: u64,
        kmercounts: Vec<KmerCount>,
        filters: &HashMap<String, String>,
    ) -> Self {
        Sketch {
            name: String::from(name),
            seq_length: Some(length),
            num_valid_kmers: Some(n_kmers),
            comment: Some(String::from("")),
            filters: Some(filters.clone()),
            hashes: kmercounts,
        }
    }

    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

    pub fn apply_filtering(&mut self, filters: &FilterParams) -> bool {
        let mut filters = filters.clone();
        let filtered_hashes = filters.filter_sketch(&self.hashes);
        let filter_stats = filters.to_serialized();
        self.hashes = filtered_hashes;
        self.filters = Some(filter_stats);
        true
    }
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
        state.serialize_field("seqLength", &self.seq_length)?;
        state.serialize_field("numValidKmers", &self.num_valid_kmers)?;
        state.serialize_field("comment", &self.comment)?;
        state.serialize_field("filters", &self.filters)?;
        state.serialize_field("hashes", &hash_list)?;
        state.serialize_field("kmers", &kmer_list)?;
        state.serialize_field("counts", &count_list)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for Sketch {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[allow(non_snake_case)]
        #[derive(Deserialize)]
        struct JSONSketch {
            pub name: String,
            pub seqLength: Option<u64>,
            pub numValidKmers: Option<u64>,
            pub comment: Option<String>,
            pub filters: Option<HashMap<String, String>>,
            hashes: Vec<QuotedU64>,
            kmers: Option<Vec<String>>,
            counts: Option<Vec<u64>>,
        }

        let mut jsketch = JSONSketch::deserialize(deserializer)?;

        let mut kmercount_list = Vec::with_capacity(jsketch.hashes.len());
        for i in 0..jsketch.hashes.len() {
            let hash = jsketch.hashes[i].0;
            let kmer = match &mut jsketch.kmers {
                Some(v) => mem::replace(&mut v[i], String::new()).into_bytes(),
                None => Vec::new(),
            };
            let count = match &jsketch.counts {
                Some(v) => v[i],
                None => 1,
            };
            kmercount_list.push(KmerCount {
                hash,
                kmer,
                count,
                extra_count: count / 2,
            });
        }
        Ok(Sketch {
            name: jsketch.name,
            seq_length: jsketch.seqLength,
            num_valid_kmers: jsketch.numValidKmers,
            comment: jsketch.comment,
            filters: jsketch.filters,
            hashes: kmercount_list,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MultiSketch {
    pub kmer: u8,
    pub alphabet: String,
    #[serde(rename = "preserveCase")]
    pub preserve_case: bool,
    pub canonical: bool,
    #[serde(rename = "sketchSize")]
    pub sketch_size: u32,
    #[serde(rename = "hashType")]
    pub hash_type: String,
    #[serde(rename = "hashBits")]
    pub hash_bits: u16,
    #[serde(rename = "hashSeed")]
    pub hash_seed: u64,
    pub scale: Option<f64>,
    pub sketches: Vec<Sketch>,
}

impl MultiSketch {
    pub fn extend(&mut self, other: &MultiSketch, name: &str) -> FinchResult<()> {
        // kmer, hashType, hashSeed, and hashBits must be same
        if self.kmer != other.kmer {
            bail!(
                "{} has a different kmer length ({}) from others ({})",
                name,
                self.kmer,
                other.kmer
            );
        } else if self.hash_type != other.hash_type {
            bail!(
                "{} used a different hash ({}) from others ({})",
                name,
                self.hash_type,
                other.hash_type
            );
        } else if self.hash_seed != other.hash_seed {
            bail!(
                "{} had a different hash seed ({}) from others ({})",
                name,
                self.hash_seed,
                other.hash_seed
            );
        } else if self.hash_bits != other.hash_bits {
            bail!(
                "{} used a different length hash ({}) from others ({})",
                name,
                self.hash_bits,
                other.hash_bits
            );
        }

        self.sketches.extend(other.sketches.clone());
        Ok(())
    }

    pub fn get_params(&self) -> FinchResult<SketchParams> {
        Ok(match (&*self.hash_type, self.scale) {
            ("MurmurHash3_x64_128", None) => {
                if self.hash_bits != 64 {
                    bail!(
                        "Multisketch has incompatible hash size ({} != 64)",
                        self.hash_bits
                    );
                }
                SketchParams::Mash {
                    kmers_to_sketch: self.sketch_size as usize,
                    final_size: self.sketch_size as usize,
                    no_strict: true,
                    kmer_length: self.kmer,
                    hash_seed: self.hash_seed,
                }
            }
            ("MurmurHash3_x64_128", Some(scale)) => {
                if self.hash_bits != 64 {
                    bail!(
                        "Multisketch has incompatible hash size ({} != 64)",
                        self.hash_bits
                    );
                }
                SketchParams::Scaled {
                    kmers_to_sketch: self.sketch_size as usize,
                    kmer_length: self.kmer,
                    scale,
                    hash_seed: self.hash_seed,
                }
            }
            ("None", _) => SketchParams::AllCounts {
                kmer_length: self.kmer,
            },
            (x, _) => bail!("{} sketch type is not supported", x),
        })
    }

    pub fn drop_all_kmers(&mut self) {
        for sketch in &mut self.sketches {
            for hash in &mut sketch.hashes {
                hash.kmer = vec![];
            }
        }
    }
}

struct QuotedU64(u64);

impl<'de> Deserialize<'de> for QuotedU64 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct QuotedU64Visitor;

        impl<'de> Visitor<'de> for QuotedU64Visitor {
            type Value = QuotedU64;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("usize as a json string")
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                value.parse().map(QuotedU64).map_err(de::Error::custom)
            }
        }

        deserializer.deserialize_str(QuotedU64Visitor)
    }
}
