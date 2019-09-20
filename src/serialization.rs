use std::collections::HashMap;
use std::convert::TryFrom;
use std::fmt;
use std::io::{BufRead, Write};
use std::mem;

#[cfg(feature = "mash_format")]
use capnp::message;
#[cfg(feature = "mash_format")]
use capnp::serialize as capnp_serialize;
use failure::{bail, format_err};
use serde::de::{self, Deserialize, Deserializer, Visitor};
use serde::ser::{Serialize, SerializeStruct, Serializer};

use crate::filtering::FilterParams;
#[cfg(feature = "mash_format")]
use crate::mash_capnp::min_hash;
use crate::sketch_schemes::{ItemHash, KmerCount, SketchParams};
use crate::Result as FinchResult;

pub const FINCH_EXT: &str = ".sk";
pub const MASH_EXT: &str = ".msh";

#[derive(Debug, Serialize, Deserialize)]
pub struct SketchDistance {
    pub containment: f64,
    pub jaccard: f64,
    #[serde(rename = "mashDistance")]
    pub mash_distance: f64,
    #[serde(rename = "commonHashes")]
    pub common_hashes: u64,
    #[serde(rename = "totalHashes")]
    pub total_hashes: u64,
    pub query: String,
    pub reference: String,
}

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
        let (filtered_hashes, low_abun) = filters.filter_sketch(&self.hashes);
        let mut filter_stats = filters.serialize_filter_params();
        if let Some(la) = low_abun {
            filter_stats.insert(String::from("minCopies"), la.to_string());
        }
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
                    scale: scale,
                    hash_seed: self.hash_seed,
                }
            }
            ("None", _) => SketchParams::AllCounts {
                kmer_length: self.kmer,
            },
            (x, _) => bail!("{} sketch type is not supported", x),
        })
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
        mash_file.set_preserve_case(sketches.preserve_case);
        mash_file.set_hash_seed(sketches.hash_seed as u32);
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
            if let Some(seq_length) = sketch.seq_length {
                mash_sketch.set_length64(seq_length);
            }
            if let Some(num_valid_kmers) = sketch.num_valid_kmers {
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
                mash_counts.reborrow().set(
                    j as u32,
                    TryFrom::try_from(hash.count).map_err(|_| {
                        format_err!("Counts greater than 32-bit not supported in mash files")
                    })?,
                );
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
        preserve_case: mash_data.get_preserve_case(),
        canonical: !mash_data.get_noncanonical(),
        sketch_size: 0,
        hash_type: String::from("MurmurHash3_x64_128"),
        hash_bits: 64u16,
        hash_seed: u64::from(mash_data.get_hash_seed()),
        scale: None,
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
                    count: u64::from(c),
                    extra_count: 0,
                })
                .collect()
        };

        sketches.sketches.push(Sketch {
            name: String::from(reference.get_name()?),
            seq_length: Some(reference.get_length64()),
            num_valid_kmers: Some(reference.get_num_valid_kmers()),
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
