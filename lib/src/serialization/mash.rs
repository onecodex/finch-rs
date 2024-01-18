use std::io::{BufRead, Write};

use capnp::message;
use capnp::serialize as capnp_serialize;

use crate::errors::FinchResult;
use crate::filtering::FilterParams;
use crate::serialization::mash_capnp::min_hash;
use crate::serialization::Sketch;
use crate::sketch_schemes::{ItemHash, KmerCount, SketchParams};

pub fn write_mash_file(file: &mut dyn Write, sketches: &[Sketch]) -> FinchResult<()> {
    let params = SketchParams::from_sketches(sketches)?;

    let mut message = message::Builder::new_default();
    {
        let mut mash_file: min_hash::Builder = message.init_root::<min_hash::Builder>();
        mash_file.set_kmer_size(u32::from(params.k()));
        mash_file.set_hash_seed(params.hash_info().2 as u32);
        mash_file.set_error(0.0); // TODO: from filters?
                                  // TODO: should we get these next 3 from a dummy method on SketchParams?
        mash_file.set_noncanonical(false);
        mash_file.set_preserve_case(false);
        mash_file.set_alphabet("ACGT");
        // not sure what these next 3 mean or if we have the right values?
        let largest_size = sketches.iter().map(|s| s.hashes.len()).max().unwrap_or(1);
        mash_file.set_window_size(u32::from(params.k()));
        mash_file.set_min_hashes_per_window(largest_size as u32);
        mash_file.set_concatenated(true);

        let mash_sketches_list = mash_file.init_reference_list();
        let mut mash_sketches = mash_sketches_list.init_references(sketches.len() as u32);

        for (i, sketch) in sketches.iter().enumerate() {
            let mut mash_sketch: min_hash::reference_list::reference::Builder =
                mash_sketches.reborrow().get(i as u32);
            mash_sketch.set_name(&sketch.name);
            mash_sketch.set_comment(&sketch.comment);
            mash_sketch.set_length64(sketch.seq_length);
            mash_sketch.set_num_valid_kmers(sketch.num_valid_kmers);
            {
                let mash_hashes = mash_sketch
                    .reborrow()
                    .init_hashes64(sketch.hashes.len() as u32);
                for (j, hash) in sketch.hashes.iter().enumerate() {
                    mash_hashes.reborrow().set(j as u32, hash.hash);
                }
            }
            let mash_counts = mash_sketch.init_counts32(sketch.hashes.len() as u32);
            for (j, hash) in sketch.hashes.iter().enumerate() {
                mash_counts.reborrow().set(j as u32, hash.count);
            }
        }
    }

    capnp_serialize::write_message(file, &message)?;
    Ok(())
}

pub fn read_mash_file(file: &mut dyn BufRead) -> FinchResult<Vec<Sketch>> {
    let options = *message::ReaderOptions::new().traversal_limit_in_words(
        // measured in words
        // 1 word = 8 bytes
        Some(2 * 1024 * 1024 * 1024),
    );
    let reader = capnp_serialize::read_message(file, options)?;
    let mash_data: min_hash::Reader = reader.get_root::<min_hash::Reader>()?;

    let kmers_to_sketch = 0;

    let sketch_params = SketchParams::Mash {
        kmers_to_sketch,
        final_size: kmers_to_sketch,
        no_strict: true,
        hash_seed: u64::from(mash_data.get_hash_seed()),
        kmer_length: mash_data.get_kmer_size() as u8,
    };

    /*
        alphabet: String::from(mash_data.get_alphabet()?),
        preserve_case: mash_data.get_preserve_case(),
        canonical: !mash_data.get_noncanonical(),
    */

    let reference_list = mash_data.get_reference_list()?;
    let reference_list_old = mash_data.get_reference_list_old()?;

    let references = if reference_list.has_references() {
        reference_list.get_references()?
    } else {
        reference_list_old.get_references()?
    };

    let mut sketches = Vec::new();
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
                    label: None,
                })
                .collect()
        } else {
            hashes
                .iter()
                .zip(counts.iter())
                .map(|(h, c)| KmerCount {
                    hash: h as ItemHash,
                    kmer: Vec::new(),
                    count: c,
                    extra_count: c / 2,
                    label: None,
                })
                .collect()
        };

        sketches.push(Sketch {
            name: String::from(reference.get_name()?),
            seq_length: reference.get_length64(),
            num_valid_kmers: reference.get_num_valid_kmers(),
            comment: String::from(reference.get_comment()?),
            hashes: kmercounts,
            sketch_params: sketch_params.clone(),
            filter_params: FilterParams::default(),
        });
    }

    Ok(sketches)
}
