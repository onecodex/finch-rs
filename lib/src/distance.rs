use std::cmp::Ordering;

#[cfg(feature = "numpy")]
use numpy::ndarray::Array2;

use crate::serialization::{Sketch, SketchDistance};
use crate::sketch_schemes::KmerCount;

pub fn distance(
    query_sketch: &Sketch,
    ref_sketch: &Sketch,
    old_mode: bool,
) -> Result<SketchDistance, &'static str> {
    let distances = if old_mode {
        old_distance(&query_sketch.hashes, &ref_sketch.hashes)
    } else {
        // since we always examine to the lowest of the sketch maxima, a
        // min_scale of 0 is a noop; otherwise we only set a scale if both of
        // the sketches are scaled (there may be a slight improvement in
        // comparing a unscaled "higher range" sketch to a scaled lower range
        // using the scale, but that makes things more complicated because we
        // need two scale values, etc)
        let mut min_scale = 0.;
        if let Some(scale1) = query_sketch.sketch_params.hash_info().3 {
            if let Some(scale2) = ref_sketch.sketch_params.hash_info().3 {
                min_scale = f64::min(scale1, scale2);
            }
        }
        raw_distance(&query_sketch.hashes, &ref_sketch.hashes, min_scale)
    };

    let containment = distances.0;
    let jaccard = distances.1;
    let common_hashes = distances.2;
    let total_hashes = distances.3;
    let k = query_sketch.sketch_params.k() as f64;
    let mash_distance: f64 = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / k;
    Ok(SketchDistance {
        containment,
        jaccard,
        mash_distance: f64::min(1f64, f64::max(0f64, mash_distance)),
        common_hashes,
        total_hashes,
        query: query_sketch.name.to_string(),
        reference: ref_sketch.name.to_string(),
    })
}

/// Estimates set statistics based on two slices of `KmerCount` sorted by hash,
/// ignoring hashes past a certain point (details below).
///
/// Returns a tuple of the form (containment,
/// [jaccard index](https://en.wikipedia.org/wiki/Jaccard_index),
/// size of [intersection](https://en.wikipedia.org/wiki/Intersection_(set_theory)),
/// size of [union](https://en.wikipedia.org/wiki/Union_(set_theory))).
///
/// Ignores hashes that are greater than `scale` * the maximum possible hash
/// (i.e. `u64::max_value()`). A scale of `0_f64` indicates to not ignore
/// hashes for this reason. Additionally ignores a slice's hashes if they are
/// greater than the other slice's maximum hash; this method better
/// approximates the document distance when the documents are of different
/// or of unknown sizes.
///
/// If the `KmerCount` slice arguments are not sorted by hash, the values of
/// the returned set statistics are unspecified.
pub fn raw_distance(
    query_hashes: &[KmerCount],
    ref_hashes: &[KmerCount],
    scale: f64,
) -> (f64, f64, u64, u64) {
    fn kmers_are_sorted(kmer_counts: &[KmerCount]) -> bool {
        for slice in kmer_counts.windows(2) {
            if slice[0].hash > slice[1].hash {
                return false;
            }
        }
        true
    }
    debug_assert!(kmers_are_sorted(query_hashes));
    debug_assert!(kmers_are_sorted(ref_hashes));

    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut common: u64 = 0;
    while let (Some(query), Some(refer)) = (query_hashes.get(i), ref_hashes.get(j)) {
        match query.hash.cmp(&refer.hash) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                common += 1;
                i += 1;
                j += 1;
            }
        }
    }

    // at this point we've exhausted one of the two sketches, but we may have
    // more counts in the other to compare if these were scaled sketches
    if scale > 0. {
        let max_hash = u64::MAX / scale.recip() as u64;
        while query_hashes
            .get(i)
            .map(|kmer_count| kmer_count.hash < max_hash)
            .unwrap_or(false)
        {
            i += 1;
        }
        while ref_hashes
            .get(j)
            .map(|kmer_count| kmer_count.hash < max_hash)
            .unwrap_or(false)
        {
            j += 1;
        }
    }

    let containment = if j == 0 { 0. } else { common as f64 / j as f64 };
    let total = i as u64 - common + j as u64;
    let jaccard: f64 = if total == 0 {
        1.
    } else {
        common as f64 / total as f64
    };

    (containment, jaccard, common, total)
}

/// This computes set statistics from one set of hashes to another.
///
/// Every hash in the reference set is considered while only those hashes in the
/// query set that are in the same range as the reference set are compared. This
/// should be a more accurate representation of the query set's containment in
/// the reference set because we consider all of the reference set. In
/// practice, there may be issues especially if the query is sketched to a
/// different effective scale than the reference.
pub fn old_distance(query_sketch: &[KmerCount], ref_sketch: &[KmerCount]) -> (f64, f64, u64, u64) {
    let mut i: usize = 0;
    let mut common: u64 = 0;
    let mut total: u64 = 0;

    for ref_hash in ref_sketch {
        while (query_sketch[i].hash < ref_hash.hash) && (i < query_sketch.len() - 1) {
            i += 1;
        }

        if query_sketch[i].hash == ref_hash.hash {
            common += 1;
        }

        total += 1;
    }

    // Numerator is A-intersect-B, |A| is the denominator, we enforce |A| == |B|
    let containment: f64 = common as f64 / total as f64;
    let jaccard: f64 = common as f64 / (common + 2 * (total - common)) as f64;
    (containment, jaccard, common, total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn kc(arr: &[u64]) -> Vec<KmerCount> {
        arr.iter()
            .map(|x| KmerCount {
                hash: *x,
                kmer: vec![],
                count: 1,
                extra_count: 1,
                label: None,
            })
            .collect()
    }

    proptest! {
        #[test]
        fn test_raw_distance_commutes(mut query_hashes: Vec<u64>, mut ref_hashes: Vec<u64>) {
            query_hashes.sort();
            ref_hashes.sort();
            let lhs = kc(&query_hashes);
            let rhs = kc(&ref_hashes);
            prop_assert_eq!(raw_distance(&lhs, &rhs, 0.), raw_distance(&rhs, &lhs, 0.));
        }
    }

    #[test]
    fn test_raw_distance() {
        let (cont, jac, com, total) = raw_distance(&kc(&[0, 1, 2]), &kc(&[1, 2]), 0.);
        assert_eq!(cont, 2. / 2.);
        assert_eq!(jac, 2. / 3.);
        assert_eq!(com, 2);
        assert_eq!(total, 3);

        let (cont, jac, com, total) = raw_distance(&kc(&[0, 2]), &kc(&[1, 2]), 0.);
        assert_eq!(cont, 1. / 2.);
        assert_eq!(jac, 1. / 3.);
        assert_eq!(com, 1);
        assert_eq!(total, 3);

        let (cont, jac, com, total) = raw_distance(&kc(&[0, 1]), &kc(&[2, 3]), 0.);
        assert_eq!(cont, 0. / 2.);
        assert_eq!(jac, 0. / 2.);
        assert_eq!(com, 0);
        assert_eq!(total, 2);

        assert_eq!((0., 1., 0, 0), raw_distance(&kc(&[]), &kc(&[]), 0.));
        assert_eq!((0., 1., 0, 0), raw_distance(&kc(&[]), &kc(&[5]), 0.));
    }

    #[test]
    fn test_raw_distance_scaled() {
        // note  below that a scale cutoff of 1e-18 translates to a max_hash of 18

        // if the hashes extend above the scale just ignore
        let (cont, jac, com, total) = raw_distance(&kc(&[10, 15, 20]), &kc(&[15, 20]), 1e-18);
        assert_eq!(cont, 2. / 2.);
        assert_eq!(jac, 2. / 3.);
        assert_eq!(com, 2);
        assert_eq!(total, 3);

        // otherwise, include the extra hashes outside the range
        let (cont, jac, com, total) = raw_distance(&kc(&[5, 10, 15]), &kc(&[5, 10]), 1e-18);
        assert_eq!(cont, 2. / 2.);
        assert_eq!(jac, 2. / 3.);
        assert_eq!(com, 2);
        assert_eq!(total, 3);

        // only include up to the scale boundary though
        let (cont, jac, com, total) = raw_distance(&kc(&[5, 10, 15, 20]), &kc(&[5, 10]), 1e-18);
        assert_eq!(cont, 2. / 2.);
        assert_eq!(jac, 2. / 3.);
        assert_eq!(com, 2);
        assert_eq!(total, 3);

        // and check in the reverse for the containment
        let (cont, jac, com, total) = raw_distance(&kc(&[5, 10]), &kc(&[5, 10, 15, 20]), 1e-18);
        assert_eq!(cont, 2. / 3.);
        assert_eq!(jac, 2. / 3.);
        assert_eq!(com, 2);
        assert_eq!(total, 3);
    }

    /// This is a straight transcription of the calculation used in Mash
    /// itself; now used for compatibility testing
    fn mash_paper_distance(sketch2: &[KmerCount], sketch1: &[KmerCount]) -> (f64, f64, u64, u64) {
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

    #[test]
    fn test_mash_compatibility() {
        // these are the same test conditions as in `test_raw_distance` but
        // against our translation of the mash distance code. this is just
        // here to sanity check some of our assumptions about how the
        // distances should be calculated; the main difference being in how
        // the original formula calculated the denominator

        let (cont, _jac, _com, _total) = mash_paper_distance(&kc(&[0, 1, 2]), &kc(&[1, 2]));
        assert_eq!(cont, 2. / 2.); // note this is actually 1./1.
                                   // assert_eq!(jac, 2. / 3.);
                                   // assert_eq!(com, 2);
                                   // assert_eq!(total, 3);

        let (_cont, _jac, _com, _total) = mash_paper_distance(&kc(&[0, 2]), &kc(&[1, 2]));
        // assert_eq!(cont, 1. / 2.);
        // assert_eq!(jac, 1. / 3.);
        // assert_eq!(com, 1);
        // assert_eq!(total, 3);

        let (_cont, jac, com, total) = mash_paper_distance(&kc(&[0, 1]), &kc(&[2, 3]));
        // assert_eq!(cont, 0. / 2.);
        assert_eq!(jac, 0. / 2.);
        assert_eq!(com, 0);
        assert_eq!(total, 2);
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
        let array1 = queue1.to_sketch();

        let mut queue2 = ScaledSketcher::new(3, 0.001, 2, 42);
        queue2.push(b"ca", 0);
        queue2.push(b"cc", 1);
        queue2.push(b"ac", 0);
        queue2.push(b"ac", 1);
        let array2 = queue2.to_sketch();

        let dist = distance(&array1, &array2, false)?;
        assert_eq!(dist.jaccard, 1.0);
        assert_eq!(dist.containment, 1.0);
        assert_eq!(dist.common_hashes, 3);

        Ok(())
    }
}

// TODO: add another method like this to allow 0's in ref sketch for hashes present in sketches?
// TODO: maybe we want to do NNLS on these matrices in Rust? for example code, see:
// https://github.com/igmanthony/fnnls/blob/master/src/fnnls.rs
// (for comments about that code also see https://github.com/rust-ndarray/ndarray/issues/649 )
#[cfg(feature = "numpy")]
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
