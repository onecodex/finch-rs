use std::collections::HashMap;
use std::fmt;
use std::io::{BufRead, Write};
use std::mem;

#[cfg(feature = "mash_format")]
use capnp::message;
#[cfg(feature = "mash_format")]
use capnp::serialize as capnp_serialize;
use serde::de::{self, Deserialize, Deserializer, Visitor};
use serde::ser::{Serialize, SerializeStruct, Serializer};

use crate::filtering::{filter_sketch, FilterParams};
#[cfg(feature = "mash_format")]
use crate::mash_capnp::min_hash;
use crate::minhashes::{ItemHash, KmerCount};
use crate::Result as FinchResult;

pub const FINCH_EXT: &str = ".sk";
pub const MASH_EXT: &str = ".msh";

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
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Sketch {
    pub name: String,
    pub seqLength: Option<u64>,
    pub numValidKmers: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    pub hashes: Vec<KmerCount>,
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
            hashes: Vec<QuotedUsize>,
            kmers: Option<Vec<String>>,
            counts: Option<Vec<u16>>,
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
    pub fn new(
        name: &str,
        length: u64,
        n_kmers: u64,
        kmercounts: Vec<KmerCount>,
        filters: &HashMap<String, String>,
    ) -> Self {
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

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

    pub fn apply_filtering(&mut self, filters: &FilterParams) -> bool {
        let (filtered_hashes, filter_stats) = filter_sketch(&self.hashes, &filters);
        self.hashes = filtered_hashes;
        self.filters = Some(filter_stats);
        true
    }
}

#[cfg(feature = "mash_format")]
pub fn write_mash_file(mut file: &mut dyn Write, sketches: &MultiSketch) -> FinchResult<()> {
    let mut message = message::Builder::new_default();
    {
        let mut mash_file: min_hash::Builder = message.init_root::<min_hash::Builder>();
        mash_file.set_kmer_size(u32::from(sketches.kmer));
        mash_file.set_window_size(u32::from(sketches.kmer));
        mash_file.set_error(0.0); // TODO: from filters?
        mash_file.set_noncanonical(!sketches.canonical);
        mash_file.set_preserve_case(sketches.preserveCase);
        mash_file.set_hash_seed(sketches.hashSeed as u32);
        mash_file.set_alphabet(&sketches.alphabet);
        // not sure what these two mean?
        let largest_size = sketches
            .sketches
            .iter()
            .map(|s| s.hashes.len())
            .max()
            .unwrap_or(1);
        mash_file.set_min_hashes_per_window(largest_size as u32);
        mash_file.set_concatenated(true);

        let mash_sketches_list = mash_file.init_reference_list();
        let mut mash_sketches = mash_sketches_list.init_references(sketches.sketches.len() as u32);

        for (i, sketch) in sketches.sketches.iter().enumerate() {
            let mut mash_sketch: min_hash::reference_list::reference::Builder =
                mash_sketches.reborrow().get(i as u32);
            mash_sketch.set_name(&sketch.name);
            if let Some(ref comment) = sketch.comment {
                mash_sketch.set_comment(&comment);
            }
            if let Some(seq_length) = sketch.seqLength {
                mash_sketch.set_length64(seq_length);
            }
            if let Some(num_valid_kmers) = sketch.numValidKmers {
                mash_sketch.set_num_valid_kmers(num_valid_kmers);
            }
            {
                let mash_hashes = mash_sketch
                    .reborrow()
                    .init_hashes64(sketch.hashes.len() as u32);
                for (j, hash) in sketch.hashes.iter().enumerate() {
                    mash_hashes.reborrow().set(j as u32, hash.hash as u64);
                }
            }
            let mash_counts = mash_sketch.init_counts32(sketch.hashes.len() as u32);
            for (j, hash) in sketch.hashes.iter().enumerate() {
                mash_counts.reborrow().set(j as u32, u32::from(hash.count));
            }
        }
    }

    capnp_serialize::write_message(&mut file, &message)?;
    Ok(())
}

#[cfg(feature = "mash_format")]
pub fn read_mash_file(mut file: &mut dyn BufRead) -> FinchResult<MultiSketch> {
    let options = *message::ReaderOptions::new().traversal_limit_in_words(64 * 1024 * 1024);
    let reader = capnp_serialize::read_message(&mut file, options)?;
    let mash_data: min_hash::Reader = reader.get_root::<min_hash::Reader>()?;

    let mut sketches = MultiSketch {
        kmer: mash_data.get_kmer_size() as u8,
        alphabet: String::from(mash_data.get_alphabet()?),
        preserveCase: mash_data.get_preserve_case(),
        canonical: !mash_data.get_noncanonical(),
        sketchSize: 0,
        hashType: String::from("MurmurHash3_x64_128"),
        hashBits: 64u16,
        hashSeed: u64::from(mash_data.get_kmer_size()),
        sketches: Vec::new(),
    };

    let reference_list = mash_data.get_reference_list()?;
    let reference_list_old = mash_data.get_reference_list_old()?;

    let references = if reference_list.has_references() {
        reference_list.get_references()?
    } else {
        reference_list_old.get_references()?
    };

    for reference in references {
        let hashes = reference.get_hashes64()?;
        let counts = reference.get_counts32()?;
        let kmercounts = if counts.len() == 0 {
            // reference_list_old doesn't seem to have counts?
            hashes
                .iter()
                .map(|h| KmerCount {
                    hash: h as ItemHash,
                    kmer: Vec::new(),
                    count: 1,
                    extra_count: 0,
                })
                .collect()
        } else {
            hashes
                .iter()
                .zip(counts.iter())
                .map(|(h, c)| KmerCount {
                    hash: h as ItemHash,
                    kmer: Vec::new(),
                    count: c as u16,
                    extra_count: 0,
                })
                .collect()
        };

        sketches.sketches.push(Sketch {
            name: String::from(reference.get_name()?),
            seqLength: Some(reference.get_length64()),
            numValidKmers: Some(reference.get_num_valid_kmers()),
            comment: Some(String::from(reference.get_comment()?)),
            filters: None,
            hashes: kmercounts,
        });
    }

    Ok(sketches)
}

#[cfg(not(feature = "mash_format"))]
pub fn write_mash_file(mut file: &mut Write, sketches: &MultiSketch) -> FinchResult<()> {
    bail!("Finch wasn't compiled with Mash format support")
}

#[cfg(not(feature = "mash_format"))]
pub fn read_mash_file(mut file: &mut BufRead) -> FinchResult<MultiSketch> {
    bail!("Finch wasn't compiled with Mash format support")
}

struct QuotedUsize(usize);

impl<'de> Deserialize<'de> for QuotedUsize {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct QuotedUsizeVisitor;

        impl<'de> Visitor<'de> for QuotedUsizeVisitor {
            type Value = QuotedUsize;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("usize as a json string")
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                value.parse().map(QuotedUsize).map_err(de::Error::custom)
            }
        }

        deserializer.deserialize_str(QuotedUsizeVisitor)
    }
}
