use std::cmp;

use crate::minhashes::KmerCount;


pub fn cardinality(sketch: &[KmerCount]) -> Result<u64, &'static str> {
    // Other (possibly more accurate) possibilities:
    // "hyper log-log" estimate from lowest value?
    // multiset distribution applied to total count number?
    // "AKMV" approach: http://people.mpi-inf.mpg.de/~rgemulla/publications/beyer07distinct.pdf

    // fast and simple k-minimum value estimate
    // https://research.neustar.biz/2012/07/09/sketch-of-the-day-k-minimum-values/
    if sketch.is_empty() {
        return Ok(0u64)
    }
    Ok(((sketch.len() - 1) as f32 / (sketch.last().unwrap().hash as f32 / usize::max_value() as f32)) as u64)
}


/// Generates a Vec of numbers of kmers for each coverage level
///
/// For example, a size 1000 sketch of the same genome repeated 5 times (e.g. 5x coverage) should
/// produce a "histogram" like [0, 0, 0, 0, 1000] (assuming no repetative kmers in the genome)
///
pub fn hist(sketch: &[KmerCount]) -> Vec<u64> {
    let mut counts = vec![0u64; 65536];
    let mut max_count: u16 = 0;
    for kmer in sketch {
        max_count = cmp::max(max_count, kmer.count);
        counts[kmer.count as usize - 1] += 1;
    }
    counts.truncate(max_count as usize);
    counts
}


#[test]
fn test_hist() {
    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 3, kmer: vec![], count: 1, extra_count: 0},
    ];

    let hist_data = hist(&sketch);
    assert_eq!(hist_data.len(), 1);
    assert_eq!(hist_data[0], 3);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 4, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 2, extra_count: 0},
        KmerCount {hash: 3, kmer: vec![], count: 4, extra_count: 0},
        KmerCount {hash: 4, kmer: vec![], count: 3, extra_count: 0},
    ];

    let hist_data = hist(&sketch);
    assert_eq!(hist_data.len(), 4);
    assert_eq!(hist_data[0], 0);
    assert_eq!(hist_data[1], 1);
    assert_eq!(hist_data[2], 1);
    assert_eq!(hist_data[3], 2);
}
