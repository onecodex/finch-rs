use ndarray::Array2;

use crate::hash_schemes::KmerCount;
use crate::serialization::SketchDistance;

pub fn distance(
    sketch1: &[KmerCount],
    sketch2: &[KmerCount],
    sketch1_name: &str,
    sketch2_name: &str,
    mash_mode: bool,
) -> Result<SketchDistance, &'static str> {
    // // TODO: in principle this is a good check, but sometimes one of the kmers will be "" if
    // // serialized without and that break this
    // if sketch1[0].kmer.len() != sketch2[0].kmer.len() {
    //     return Err("Sketches have different sized kmers");
    // }
    let distances = if mash_mode {
        raw_mash_distance(sketch1, sketch2)
    } else {
        raw_distance(sketch1, sketch2)
    };
    let containment = distances.0;
    let jaccard = distances.1;
    let common = distances.2;
    let total = distances.3;
    let mash_distance: f64 =
        -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / sketch1[0].kmer.len() as f64;
    Ok(SketchDistance {
        containment,
        jaccard,
        mashDistance: f64::min(1f64, f64::max(0f64, mash_distance)),
        commonHashes: common,
        totalHashes: total,
        query: sketch1_name.to_string(),
        reference: sketch2_name.to_string(),
    })
}

fn raw_mash_distance(sketch1: &[KmerCount], sketch2: &[KmerCount]) -> (f64, f64, u64, u64) {
    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut common: u64 = 0;
    let mut total: u64 = 0;
    let sketch_size = sketch1.len();

    while (total < sketch_size as u64) && (i < sketch1.len()) && (j < sketch2.len()) {
        if sketch1[i].hash < sketch2[j].hash {
            i += 1;
        } else if sketch2[j].hash < sketch1[i].hash {
            j += 1;
        } else {
            i += 1;
            j += 1;
            common += 1;
        }
        total += 1;
    }

    if total < sketch_size as u64 {
        if i < sketch1.len() {
            total += (sketch1.len() - 1) as u64;
        }

        if j < sketch2.len() {
            total += (sketch2.len() - 1) as u64;
        }

        if total > sketch_size as u64 {
            total = sketch_size as u64;
        }
    }

    let containment: f64 = common as f64 / i as f64;
    let jaccard: f64 = common as f64 / total as f64;
    (containment, jaccard, common, total)
}

pub fn raw_distance(sketch1: &[KmerCount], sketch2: &[KmerCount]) -> (f64, f64, u64, u64) {
    let mut j: usize = 0;
    let mut common: u64 = 0;
    let mut total: u64 = 0;

    for hash1 in sketch1 {
        while (sketch2[j].hash < hash1.hash) && (j < sketch2.len() - 1) {
            j += 1;
        }

        if sketch2[j].hash == hash1.hash {
            common += 1;
        }

        total += 1;
    }

    // Numerator is A-intersect-B, |A| is the denominator, we enforce |A| == |B|
    let containment: f64 = common as f64 / total as f64;
    let jaccard: f64 = common as f64 / (common + 2 * (total - common)) as f64;
    (containment, jaccard, common, total)
}

pub fn common_counts(sketch1: &[KmerCount], sketch2: &[KmerCount]) -> (u64, u64, u64, u64, u64) {
    let mut common: u64 = 0;
    let mut pos1: usize = 0;
    let mut pos2: usize = 0;
    let mut count1: u64 = 0;
    let mut count2: u64 = 0;

    while (pos1 < sketch1.len()) && (pos2 < sketch2.len()) {
        if sketch1[pos1].hash < sketch2[pos2].hash {
            pos1 += 1;
        } else if sketch2[pos2].hash < sketch1[pos1].hash {
            pos2 += 1;
        } else {
            count1 += sketch1[pos1].count;
            count2 += sketch2[pos2].count;
            pos1 += 1;
            pos2 += 1;
            common += 1;
        }
    }

    (common, pos1 as u64, pos2 as u64, count1, count2)
}

// TODO: add another method like this to allow 0's in ref sketch for hashes present in sketches?
pub fn minmer_matrix<U>(ref_sketch: &[KmerCount], sketches: &[U]) -> Array2<u64>
where
    U: AsRef<[KmerCount]>,
{
    let mut result = Array2::<u64>::zeros((sketches.len(), ref_sketch.len()));

    for (i, sketch) in sketches.iter().map(|s| s.as_ref()).enumerate() {
        let mut ref_pos = 0;
        for hash in sketch.iter() {
            while (hash.hash > ref_sketch[ref_pos].hash) && (ref_pos < ref_sketch.len() - 1) {
                ref_pos += 1;
            }

            if hash.hash == ref_sketch[ref_pos].hash {
                result[[i, ref_pos]] = hash.count;
            }
        }
    }
    result
}
