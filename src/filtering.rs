use std::cmp;
use std::collections::HashMap;

use minhashes::KmerCount;
use statistics::hist;


/// Used to pass around filter options for sketching
#[derive(Debug)]
pub struct FilterParams {
    pub filter_on: Option<bool>,
    pub abun_filter: (Option<u16>, Option<u16>),
    pub err_filter: f32,
    pub strand_filter: f32,
}

/// Applies filter options to a sketch
pub fn filter_sketch(hashes: &[KmerCount], filters: &FilterParams) -> (Vec<KmerCount>, HashMap<String, String>) {
    let filter_on = filters.filter_on.expect("Sorry! Filter should have either been passed or set during detection");
    let mut filter_stats: HashMap<String, String> = HashMap::new();

    let mut low_abun_filter = filters.abun_filter.0;

    if filter_on && filters.err_filter > 0f32 {
        let cutoff = guess_filter_threshold(&hashes, filters.err_filter);
        if let None = low_abun_filter {
            low_abun_filter = Some(cutoff);
            filter_stats.insert(String::from("errFilter"), filters.err_filter.to_string());
        }
    }

    let mut filtered_hashes = hashes.to_vec();
    if filter_on && (low_abun_filter != None || filters.abun_filter.1 != None) {
        filtered_hashes = filter_abundance(&filtered_hashes, low_abun_filter, filters.abun_filter.1);
        if let Some(v) = low_abun_filter {
            filter_stats.insert(String::from("minCopies"), v.to_string());
        }
        if let Some(v) = filters.abun_filter.1 {
            filter_stats.insert(String::from("maxCopies"), v.to_string());
        }
    }

    if filter_on && filters.strand_filter > 0f32 {
        filtered_hashes = filter_strands(&filtered_hashes, filters.strand_filter);
        filter_stats.insert(String::from("strandFilter"), filters.strand_filter.to_string());
    }

    (filtered_hashes, filter_stats)
}


/// Determines a dynamic filtering threshold for low abundance kmers
///
/// Useful for removing, e.g. kmers containing sequencing errors
///
pub fn guess_filter_threshold(sketch: &[KmerCount], filter_level: f32) -> u16 {
    let hist_data = hist(sketch);
    let total_counts = hist_data.iter().enumerate().map(|t| (t.0 as u64 + 1) * t.1).sum::<u64>() as f32;
    let cutoff_amt = filter_level * total_counts;

    // calculate the coverage that N% of the weighted data is above
    let mut wgt_cutoff: usize = 1;
    let mut cum_count: u64 = 0;
    for count in &hist_data {
        cum_count += wgt_cutoff as u64 * *count as u64;
        if cum_count as f32 > cutoff_amt {
            break;
        }
        wgt_cutoff += 1;
    }

    if wgt_cutoff <= 2 {
        return wgt_cutoff as u16;
    }

    // now find the right-most global maxima
    let win_size = cmp::max(1, wgt_cutoff / 20);
    let mut sum: u64 = hist_data[..win_size].iter().sum();
    let mut lowest_val = sum;
    let mut lowest_idx = win_size;
    for (i, j) in (0..wgt_cutoff - win_size).zip(win_size..wgt_cutoff) {
        if sum <= lowest_val {
            lowest_val = sum;
            lowest_idx = j;
        }
        sum -= hist_data[i];
        sum += hist_data[j];
    }

    lowest_idx as u16
}

#[test]
fn test_guess_filter_threshold() {
    let sketch = vec![];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 1);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 1);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 1, extra_count: 0},
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 1);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 9, extra_count: 0},
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 8);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 10, extra_count: 0},
        KmerCount {hash: 3, kmer: vec![], count: 10, extra_count: 0},
        KmerCount {hash: 4, kmer: vec![], count: 9, extra_count: 0},
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.1);
    assert_eq!(cutoff, 8);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 3, kmer: vec![], count: 2, extra_count: 0},
        KmerCount {hash: 4, kmer: vec![], count: 4, extra_count: 0},
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.1);
    assert_eq!(cutoff, 1);
}


pub fn filter_abundance(sketch: &[KmerCount], low: Option<u16>, high: Option<u16>) -> Vec<KmerCount> {
    let mut filtered = Vec::new();
    let lo_threshold = low.unwrap_or(0u16);
    let hi_threshold = high.unwrap_or(u16::max_value());
    for kmer in sketch {
        if lo_threshold <= kmer.count && kmer.count <= hi_threshold {
            filtered.push(kmer.clone());
        }
    }
    filtered
}

#[test]
fn test_filter_abundance() {
    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 1, extra_count: 0},
    ];
    let filtered = filter_abundance(&sketch, Some(1), None);
    assert_eq!(filtered.len(), 2);
    assert_eq!(filtered[0].hash, 1);
    assert_eq!(filtered[1].hash, 2);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 1, extra_count: 0},
        KmerCount {hash: 2, kmer: vec![], count: 10, extra_count: 0},
        KmerCount {hash: 3, kmer: vec![], count: 10, extra_count: 0},
        KmerCount {hash: 4, kmer: vec![], count: 9, extra_count: 0},
    ];
    let filtered = filter_abundance(&sketch, Some(9), None);
    assert_eq!(filtered.len(), 3);
    assert_eq!(filtered[0].hash, 2);
    assert_eq!(filtered[1].hash, 3);
    assert_eq!(filtered[2].hash, 4);

    let filtered = filter_abundance(&sketch, Some(2), Some(9));
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].hash, 4);
}


/// Filter out kmers that have a large abundance difference between being seen in the
/// "forward" and "reverse" orientations (picked arbitrarily which is which).
///
/// These tend to be sequencing adapters.
pub fn filter_strands(sketch: &[KmerCount], ratio_cutoff: f32) -> Vec<KmerCount> {
    let mut filtered = Vec::new();
    for kmer in sketch {
        // "special-case" anything with fewer than 16 kmers -> these are too stochastic to accurately
        // determine if something is an adapter or not. The odds of randomly picking less than 10%
        // (0 or 1 reversed kmers) in 16 should be ~ 17 / 2 ** 16 or 1/4000 so we're avoiding
        // removing "good" kmers
        if kmer.count < 16u16 {
            filtered.push(kmer.clone());
            continue;
        }

        // check the forward/reverse ratio and only add if it's within bounds
        let lowest_strand_count: u16 = cmp::min(kmer.extra_count, kmer.count - kmer.extra_count);
        if (lowest_strand_count as f32 / kmer.count as f32) >= ratio_cutoff {
            filtered.push(kmer.clone());
        }
    }
    filtered
}


#[test]
fn test_filter_strands() {
    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 10, extra_count: 1},
        KmerCount {hash: 2, kmer: vec![], count: 10, extra_count: 2},
        KmerCount {hash: 3, kmer: vec![], count: 10, extra_count: 8},
        KmerCount {hash: 4, kmer: vec![], count: 10, extra_count: 9},
    ];
    let filtered = filter_strands(&sketch, 0.15);
    assert_eq!(filtered.len(), 4);
    assert_eq!(filtered[0].hash, 1);
    assert_eq!(filtered[3].hash, 4);

    let sketch = vec![
        KmerCount {hash: 1, kmer: vec![], count: 16, extra_count: 1},
        KmerCount {hash: 2, kmer: vec![], count: 16, extra_count: 2},
        KmerCount {hash: 3, kmer: vec![], count: 16, extra_count: 8},
        KmerCount {hash: 4, kmer: vec![], count: 16, extra_count: 9},
    ];
    let filtered = filter_strands(&sketch, 0.15);
    assert_eq!(filtered.len(), 2);
    assert_eq!(filtered[0].hash, 3);
    assert_eq!(filtered[1].hash, 4);
}
