use std::convert::TryFrom;
use std::io::{BufRead, Write};

use capnp::message;
use capnp::serialize as capnp_serialize;
use failure::format_err;

use crate::serialization::mash_capnp::min_hash;
use crate::serialization::{MultiSketch, Sketch};
use crate::sketch_schemes::{ItemHash, KmerCount};
use crate::Result as FinchResult;

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
