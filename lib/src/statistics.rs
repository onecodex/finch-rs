use std::cmp;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;

use crate::sketch_schemes::hashing::NoHashHasher;
use crate::sketch_schemes::KmerCount;

pub fn cardinality(sketch: &[KmerCount]) -> Result<u64, &'static str> {
    // Other (possibly more accurate) possibilities:
    // "hyper log-log" estimate from lowest value?
    // multiset distribution applied to total count number?
    // "AKMV" approach: http://people.mpi-inf.mpg.de/~rgemulla/publications/beyer07distinct.pdf

    // fast and simple k-minimum value estimate
    // https://research.neustar.biz/2012/07/09/sketch-of-the-day-k-minimum-values/
    if sketch.is_empty() {
        return Ok(0u64);
    }
    Ok(
        ((sketch.len() - 1) as f32 / (sketch.last().unwrap().hash as f32 / usize::MAX as f32))
            as u64,
    )
}

/// Generates a Vec of numbers of kmers for each coverage level
///
/// For example, a size 1000 sketch of the same genome repeated 5 times (e.g. 5x coverage) should
/// produce a "histogram" like [0, 0, 0, 0, 1000] (assuming no repetitive kmers in the genome)
///
pub fn hist(sketch: &[KmerCount]) -> Vec<u64> {
    let mut max_count = 0;
    let mut counts: HashMap<usize, u64, BuildHasherDefault<NoHashHasher>> =
        HashMap::with_capacity_and_hasher(150_000, BuildHasherDefault::default());
    for kmer in sketch {
        max_count = cmp::max(max_count, u64::from(kmer.count));
        counts
            .entry(kmer.count as usize - 1)
            .and_modify(|c| *c += 1)
            .or_insert(1);
    }

    let mut counts_vec: Vec<u64> = Vec::with_capacity(max_count as usize);
    for i in 0..max_count {
        counts_vec.push(*counts.get(&(i as usize)).unwrap_or(&0))
    }
    counts_vec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hist() {
        let sketch = vec![
            KmerCount {
                hash: 1,
                kmer: vec![],
                count: 1,
                extra_count: 0,
                label: None,
            },
            KmerCount {
                hash: 2,
                kmer: vec![],
                count: 1,
                extra_count: 0,
                label: None,
            },
            KmerCount {
                hash: 3,
                kmer: vec![],
                count: 1,
                extra_count: 0,
                label: None,
            },
        ];

        let hist_data = hist(&sketch);
        assert_eq!(hist_data.len(), 1);
        assert_eq!(hist_data[0], 3);

        let sketch = vec![
            KmerCount {
                hash: 1,
                kmer: vec![],
                count: 4,
                extra_count: 0,
                label: None,
            },
            KmerCount {
                hash: 2,
                kmer: vec![],
                count: 2,
                extra_count: 0,
                label: None,
            },
            KmerCount {
                hash: 3,
                kmer: vec![],
                count: 4,
                extra_count: 0,
                label: None,
            },
            KmerCount {
                hash: 4,
                kmer: vec![],
                count: 3,
                extra_count: 0,
                label: None,
            },
            // https://github.com/onecodex/finch-rs/issues/63
            KmerCount {
                hash: 3,
                kmer: vec![],
                count: 126497,
                extra_count: 0,
                label: None,
            },
        ];

        let hist_data = hist(&sketch);
        assert_eq!(hist_data.len(), 126497);
        assert_eq!(hist_data[0], 0);
        assert_eq!(hist_data[1], 1);
        assert_eq!(hist_data[2], 1);
        assert_eq!(hist_data[3], 2);
        assert_eq!(hist_data[126497 - 1], 1);
    }
}
