use ndarray::Array2;

use crate::serialization::SketchDistance;
use crate::sketch_schemes::KmerCount;

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
        mash_distance: f64::min(1f64, f64::max(0f64, mash_distance)),
        common_hashes: common,
        total_hashes: total,
        query: sketch1_name.to_string(),
        reference: sketch2_name.to_string(),
    })
}

pub fn distance_scaled(
    sketch1: &[KmerCount],
    sketch2: &[KmerCount],
    sketch1_name: &str,
    sketch2_name: &str,
) -> Result<SketchDistance, &'static str> {
    // // TODO: in principle this is a good check, but sometimes one of the kmers will be "" if
    // // serialized without and that break this
    // if sketch1[0].kmer.len() != sketch2[0].kmer.len() {
    //     return Err("Sketches have different sized kmers");
    // }

    // TODO: check if both sketches are in the same scaled factor
    let distances = raw_mash_distance(sketch1, sketch2);

    let containment = distances.0;
    let jaccard = distances.1;
    let common = distances.2;
    let total = distances.3;
    let mash_distance: f64 =
        -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / sketch1[0].kmer.len() as f64;
    Ok(SketchDistance {
        containment,
        jaccard,
        mash_distance: f64::min(1f64, f64::max(0f64, mash_distance)),
        common_hashes: common,
        total_hashes: total,
        query: sketch1_name.to_string(),
        reference: sketch2_name.to_string(),
    })
}

/// This computes the set statistics between two sets of hashes.
///
/// It stops once one of the sets has "run out" of hashes, i.e. at the
/// smallest max hash of the two sets. In general this is a better
/// approximation of the true document distance when either the two
/// original documents were of different sizes or when the two documents
/// were hashed in an unscaled fashion.
fn raw_mash_distance(sketch1: &[KmerCount], sketch2: &[KmerCount]) -> (f64, f64, u64, u64) {
    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut common: u64 = 0;
    let mut total: u64 = 0;
    let sketch_size = sketch1.len(); // is this true?

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

// FIXME!!!!
pub fn best_distance(
    sketch1: &[KmerCount],
    sketch2: &[KmerCount],
    max_min: Option<u64>,
) -> (f64, f64, u64, u64) {
    // TODO: should we be normalizing the two count sets to have equivalent range?
    // it seems like our methid here is probably not scale-robust
    // should we normalize to i/j? or s1_* and s2_* to have the same sum?
    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut s1_intersect: u64 = 0;
    let mut s1_complement: u64 = 0;
    let mut s2_intersect: u64 = 0;
    let mut s2_complement: u64 = 0;
    let mut n_intersect: u64 = 0;

    // the maximum allowable hash value; either from the one passed in or from
    // the "shortest" of the two sketches if none was passed in or if the one
    // passed in was too big
    // FIXME: if one of these arrays is empty, this will panic
    let max_range = u64::min(sketch1.last().unwrap().hash, sketch2.last().unwrap().hash);
    let max_min = u64::min(max_min.unwrap_or_else(|| max_range), max_range);

    while (sketch1[i].hash <= max_min) && (sketch2[j].hash <= max_min) {
        if sketch1[i].hash < sketch2[j].hash {
            s1_complement += u64::from(sketch1[i].count);
            i += 1;
        } else if sketch2[j].hash < sketch1[i].hash {
            s2_complement += u64::from(sketch2[j].count);
            j += 1;
        } else {
            // we could update both total and common here, but since the sum
            // is the same we just push that to the end
            s1_intersect += u64::from(sketch1[i].count);
            s2_intersect += u64::from(sketch2[j].count);
            n_intersect += 1;
            i += 1;
            j += 1;
        }
    }
    let containment: f64 = s1_intersect as f64 / (s1_intersect as f64 + s1_complement as f64);

    let common = (s1_intersect + s2_intersect) / 2;
    let total = s1_complement + common + s2_complement;
    let jaccard: f64 = common as f64 / total as f64;
    (containment, jaccard, common, total)
    // Ioffe @ Google has a paper on weighted MinHash that might be useful to understand?
    // http://static.googleusercontent.com/media/research.google.com/en/us/pubs/archive/36928.pdf
}

/// This computes set statistics from one set of hashes to another.
///
/// Every hash in the first set is considered while only those hashes in the
/// second set that are in the same range as the first set are compared. This
/// should be a more accurate representation of the second set's containment in
/// the first set because we consider all of the first set.
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

// TODO: add another method like this to allow 0's in ref sketch for hashes present in sketches?
// TODO: maybe we want to do NNLS on these matrices here too? https://github.com/igmanthony/fnnls/blob/master/src/fnnls.rs
// (for comments about that code also see https://github.com/rust-ndarray/ndarray/issues/649 )
pub fn minmer_matrix<U>(ref_sketch: &[KmerCount], sketches: &[U]) -> Array2<i32>
where
    U: AsRef<[KmerCount]>,
{
    let mut result = Array2::<i32>::zeros((sketches.len(), ref_sketch.len()));

    for (i, sketch) in sketches.iter().map(|s| s.as_ref()).enumerate() {
        let mut ref_pos = 0;
        for hash in sketch.iter() {
            while (hash.hash > ref_sketch[ref_pos].hash) && (ref_pos < ref_sketch.len() - 1) {
                ref_pos += 1;
            }

            if hash.hash == ref_sketch[ref_pos].hash {
                result[[i, ref_pos]] = hash.count as i32;
            }
        }
    }
    result
}

#[test]
fn test_distance_scaled() -> Result<(), Box<dyn std::error::Error>> {
    use crate::sketch_schemes::scaled::ScaledSketcher;
    use crate::sketch_schemes::SketchScheme;

    let mut queue1 = ScaledSketcher::new(3, 0.001, 2, 42);
    queue1.push(b"ca", 0);
    queue1.push(b"cc", 1);
    queue1.push(b"ac", 0);
    queue1.push(b"ac", 1);
    let array1 = queue1.to_vec();

    let mut queue2 = ScaledSketcher::new(3, 0.001, 2, 42);
    queue2.push(b"ca", 0);
    queue2.push(b"cc", 1);
    queue2.push(b"ac", 0);
    queue2.push(b"ac", 1);
    let array2 = queue2.to_vec();

    let dist = distance_scaled(&array1, &array2, "", "")?;
    assert_eq!(dist.jaccard, 1.0);
    assert_eq!(dist.containment, 1.0);
    assert_eq!(dist.common_hashes, 3);

    Ok(())
}
