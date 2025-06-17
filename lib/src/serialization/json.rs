use std::collections::HashMap;
use std::fmt;
use std::mem;

use serde::de::{self, Deserializer, Visitor};
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};

use crate::bail;
use crate::errors::FinchResult;
use crate::filtering::FilterParams;
use crate::serialization::Sketch;
use crate::sketch_schemes::{KmerCount, SketchParams};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct JsonSketch {
    pub name: String,
    pub seq_length: Option<u64>,
    pub num_valid_kmers: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    pub hashes: Vec<KmerCount>,
}

impl JsonSketch {
    pub fn new(
        name: &str,
        length: u64,
        n_kmers: u64,
        kmercounts: Vec<KmerCount>,
        filters: &HashMap<String, String>,
    ) -> Self {
        JsonSketch {
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
}

impl From<Sketch> for JsonSketch {
    fn from(s: Sketch) -> Self {
        JsonSketch::new(
            &s.name,
            s.seq_length,
            s.num_valid_kmers,
            s.hashes,
            &s.filter_params.to_serialized(),
        )
    }
}

impl Serialize for JsonSketch {
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

impl<'de> Deserialize<'de> for JsonSketch {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[allow(non_snake_case)]
        #[derive(Deserialize)]
        struct BaseJsonSketch {
            pub name: String,
            pub seqLength: Option<u64>,
            pub numValidKmers: Option<u64>,
            pub comment: Option<String>,
            pub filters: Option<HashMap<String, String>>,
            hashes: Vec<QuotedU64>,
            kmers: Option<Vec<String>>,
            counts: Option<Vec<u32>>,
        }

        let mut jsketch = BaseJsonSketch::deserialize(deserializer)?;

        let mut kmercount_list = Vec::with_capacity(jsketch.hashes.len());
        for i in 0..jsketch.hashes.len() {
            let hash = jsketch.hashes[i].0;
            let kmer = match &mut jsketch.kmers {
                Some(v) => mem::take(&mut v[i]).into_bytes(),
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
                label: None,
            });
        }
        Ok(JsonSketch {
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
    pub sketches: Vec<JsonSketch>,
}

impl MultiSketch {
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

    pub fn from_sketches(sketches: &[Sketch]) -> FinchResult<Self> {
        let json_sketches: Vec<JsonSketch> = sketches.iter().map(|x| (*x).clone().into()).collect();
        let sketch_params = SketchParams::from_sketches(sketches)?;
        // TODO: the scale isn't actually harmonized between the sketches at
        // this point; it probably should be?
        let (hash_type, hash_bits, hash_seed, scale) = sketch_params.hash_info();
        Ok(MultiSketch {
            alphabet: "ACGT".to_string(),
            preserve_case: false,
            canonical: true,

            sketch_size: sketch_params.expected_size() as u32,
            kmer: sketch_params.k(),
            hash_type: hash_type.to_string(),
            hash_bits,
            hash_seed,
            scale,
            sketches: json_sketches,
        })
    }

    pub fn to_sketches(&self) -> FinchResult<Vec<Sketch>> {
        let empty_hashmap = HashMap::new();
        let mut sketches = Vec::with_capacity(self.sketches.len());
        let sketch_params = self.get_params()?;
        for sketch in &self.sketches {
            let filters = sketch.filters.as_ref().unwrap_or(&empty_hashmap);
            let filter_params = FilterParams::from_serialized(filters)?;
            sketches.push(Sketch {
                name: sketch.name.clone(),
                seq_length: sketch.seq_length.unwrap_or(0),
                num_valid_kmers: sketch.num_valid_kmers.unwrap_or(0),
                comment: sketch.comment.clone().unwrap_or_default(),
                hashes: sketch.hashes.clone(),
                filter_params,
                sketch_params: sketch_params.clone(),
            });
        }
        Ok(sketches)
    }
}

struct QuotedU64(u64);

impl<'de> Deserialize<'de> for QuotedU64 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct QuotedU64Visitor;

        impl Visitor<'_> for QuotedU64Visitor {
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
