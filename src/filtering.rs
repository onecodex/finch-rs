use std::cmp;

use minhashes::{KmerCount};

/// Filters low abundance kmers out of a sketch
///
/// Useful for removing, e.g. kmers containing sequencing errors
///
pub fn filter_sketch(sketch: &[KmerCount], filter_level: f32) -> (Vec<KmerCount>, u16) {
    let hist_data = hist(sketch);
    let total_counts = hist_data.iter().enumerate().map(|t| (t.0 as u64 + 1) * t.1).sum::<u64>() as f32;
    let cutoff_amt = filter_level * total_counts;

    // calculate the first "counts" bin to keep
    let mut cutoff: u16 = 1;
    let mut cum_count: u64 = 0;
    for count in hist_data {
        cum_count += cutoff as u64 * count as u64;
        if cum_count as f32 > cutoff_amt {
            break;
        }
        cutoff += 1;
    }

    filter_abs(sketch, cutoff)
}

pub fn filter_abs(sketch: &[KmerCount], cov_cutoff: u16) -> (Vec<KmerCount>, u16) {
    let mut filtered = Vec::new();
    for kmer in sketch {
        if kmer.count >= cov_cutoff {
            filtered.push(kmer.clone());
        }
    }

    (filtered, cov_cutoff)
}


#[test]
fn test_filter_sketch() {
    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1},
        KmerCount {hash: 2, kmer: vec![], count: 1},
    ];
    let (filtered, cutoff) = filter_sketch(&sketch, 0.2);
    assert_eq!(filtered.len(), 2);
    assert_eq!(cutoff, 1);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1},
        KmerCount {hash: 2, kmer: vec![], count: 9},
    ];
    let (filtered, cutoff) = filter_sketch(&sketch, 0.2);
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].hash, 2);
    assert_eq!(cutoff, 9);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1},
        KmerCount {hash: 2, kmer: vec![], count: 10},
        KmerCount {hash: 3, kmer: vec![], count: 10},
        KmerCount {hash: 4, kmer: vec![], count: 9},
    ];
    let (filtered, cutoff) = filter_sketch(&sketch, 0.1);
    assert_eq!(filtered.len(), 3);
    assert_eq!(filtered[0].hash, 2);
    assert_eq!(filtered[1].hash, 3);
    assert_eq!(cutoff, 9);
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
        KmerCount {hash: 1, kmer: vec![], count: 1},
        KmerCount {hash: 2, kmer: vec![], count: 1},
        KmerCount {hash: 3, kmer: vec![], count: 1},
    ];

    let hist_data = hist(&sketch);
    assert_eq!(hist_data.len(), 1);
    assert_eq!(hist_data[0], 3);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 4},
        KmerCount {hash: 2, kmer: vec![], count: 2},
        KmerCount {hash: 3, kmer: vec![], count: 4},
        KmerCount {hash: 4, kmer: vec![], count: 3},
    ];

    let hist_data = hist(&sketch);
    assert_eq!(hist_data.len(), 4);
    assert_eq!(hist_data[0], 0);
    assert_eq!(hist_data[1], 1);
    assert_eq!(hist_data[2], 1);
    assert_eq!(hist_data[3], 2);
}
