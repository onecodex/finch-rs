use std::cmp;
use std::collections::HashMap;

use crate::sketch_schemes::KmerCount;
use crate::statistics::hist;

/// Used to pass around filter options for sketching
#[derive(Clone, Debug)]
pub struct FilterParams {
    pub filter_on: Option<bool>,
    pub abun_filter: (Option<u64>, Option<u64>),
    pub err_filter: f64,
    pub strand_filter: f64,
}

impl FilterParams {
    /// Returns the filtered kmer counts and the low abundance cutoff, if there
    /// was a different one determined from the err_filter.
    pub fn filter_sketch(&mut self, hashes: &[KmerCount]) -> Vec<KmerCount> {
        let filter_on = self.filter_on == Some(true);
        let mut filtered_hashes = hashes.to_vec();

        if filter_on && self.strand_filter > 0f64 {
            filtered_hashes = filter_strands(&filtered_hashes, self.strand_filter);
        }

        if filter_on && self.err_filter > 0f64 {
            let cutoff = guess_filter_threshold(&filtered_hashes, self.err_filter);
            if let Some(v) = self.abun_filter.0 {
                // there's an existing filter so we only use this one if it's stricter
                if cutoff > v {
                    self.abun_filter.0 = Some(cutoff);
                }
            } else {
                // no filter set so just use the one we determined
                self.abun_filter.0 = Some(cutoff);
            }
        }

        if filter_on && (self.abun_filter.0.is_some() || self.abun_filter.1.is_some()) {
            filtered_hashes =
                filter_abundance(&filtered_hashes, self.abun_filter.0, self.abun_filter.1);
        }

        filtered_hashes
    }

    pub fn to_serialized(&self) -> HashMap<String, String> {
        let mut filter_stats: HashMap<String, String> = HashMap::new();
        if self.filter_on != Some(true) {
            return filter_stats;
        }

        if self.strand_filter > 0f64 {
            filter_stats.insert(String::from("strandFilter"), self.strand_filter.to_string());
        }
        if self.err_filter > 0f64 {
            filter_stats.insert(String::from("errFilter"), self.err_filter.to_string());
        }
        if let Some(v) = self.abun_filter.0 {
            filter_stats.insert(String::from("minCopies"), v.to_string());
        }
        if let Some(v) = self.abun_filter.1 {
            filter_stats.insert(String::from("maxCopies"), v.to_string());
        }
        filter_stats
    }

    pub fn from_serialized(filters: &HashMap<String, String>) -> Self {
        // TODO: remove unwraps and make this return a result
        FilterParams {
            filter_on: Some(!filters.is_empty()),
            abun_filter: (
                filters.get("minCopies").map(|x| x.parse().unwrap()),
                filters.get("maxCopies").map(|x| x.parse().unwrap()),
            ),
            err_filter: filters
                .get("errFilter")
                .unwrap_or(&"0".to_string())
                .parse()
                .unwrap(),
            strand_filter: filters
                .get("strandFilter")
                .unwrap_or(&"0".to_string())
                .parse()
                .unwrap(),
        }
    }
}

impl Default for FilterParams {
    fn default() -> Self {
        FilterParams {
            filter_on: Some(false),
            abun_filter: (None, None),
            err_filter: 0.,
            strand_filter: 0.,
        }
    }
}

/// Determines a dynamic filtering threshold for low abundance kmers. The
/// cutoff returned is the lowest number of counts that should be included
/// in any final results.
///
/// Useful for removing, e.g. low-abundance kmers arising from sequencing
/// errors
///
pub fn guess_filter_threshold(sketch: &[KmerCount], filter_level: f64) -> u64 {
    let hist_data = hist(sketch);
    let total_counts = hist_data
        .iter()
        .enumerate()
        .map(|t| (t.0 as u64 + 1) * t.1)
        .sum::<u64>() as f64;
    let cutoff_amt = filter_level * total_counts;

    // calculate the coverage that N% of the weighted data is above
    // note wgt_cutoff is an index now *not* a number of counts
    let mut wgt_cutoff: usize = 0;
    let mut cum_count: u64 = 0;
    for count in &hist_data {
        cum_count += wgt_cutoff as u64 * *count as u64;
        if cum_count as f64 > cutoff_amt {
            break;
        }
        wgt_cutoff += 1;
    }

    // special case if the cutoff is the first value
    if wgt_cutoff == 0 {
        return 1;
    }

    // now find the minima within the window to the left
    let win_size = cmp::max(1, wgt_cutoff / 20);
    let mut sum: u64 = hist_data[..win_size].iter().sum();
    let mut lowest_val = sum;
    let mut lowest_idx = win_size - 1;
    for (i, j) in (0..wgt_cutoff - win_size).zip(win_size..wgt_cutoff) {
        if sum <= lowest_val {
            lowest_val = sum;
            lowest_idx = j;
        }
        sum -= hist_data[i];
        sum += hist_data[j];
    }

    lowest_idx as u64 + 1
}

#[test]
fn test_guess_filter_threshold() {
    let sketch = vec![];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 1);

    let sketch = vec![KmerCount {
        hash: 1,
        kmer: vec![],
        count: 1,
        extra_count: 0,
    }];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 1);

    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 1);

    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 9,
            extra_count: 0,
        },
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.2);
    assert_eq!(cutoff, 8);

    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 10,
            extra_count: 0,
        },
        KmerCount {
            hash: 3,
            kmer: vec![],
            count: 10,
            extra_count: 0,
        },
        KmerCount {
            hash: 4,
            kmer: vec![],
            count: 9,
            extra_count: 0,
        },
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.1);
    assert_eq!(cutoff, 8);

    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 3,
            kmer: vec![],
            count: 2,
            extra_count: 0,
        },
        KmerCount {
            hash: 4,
            kmer: vec![],
            count: 4,
            extra_count: 0,
        },
    ];
    let cutoff = guess_filter_threshold(&sketch, 0.1);
    assert_eq!(cutoff, 1);

    // check that we don't overflow
    let sketch = vec![KmerCount {
        hash: 2,
        kmer: vec![],
        count: 2,
        extra_count: 0,
    }];
    let cutoff = guess_filter_threshold(&sketch, 1.);
    assert_eq!(cutoff, 2);
}

pub fn filter_abundance(
    sketch: &[KmerCount],
    low: Option<u64>,
    high: Option<u64>,
) -> Vec<KmerCount> {
    let mut filtered = Vec::new();
    let lo_threshold = low.unwrap_or(0u64);
    let hi_threshold = high.unwrap_or(u64::max_value());
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
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
    ];
    let filtered = filter_abundance(&sketch, Some(1), None);
    assert_eq!(filtered.len(), 2);
    assert_eq!(filtered[0].hash, 1);
    assert_eq!(filtered[1].hash, 2);

    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 1,
            extra_count: 0,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 10,
            extra_count: 0,
        },
        KmerCount {
            hash: 3,
            kmer: vec![],
            count: 10,
            extra_count: 0,
        },
        KmerCount {
            hash: 4,
            kmer: vec![],
            count: 9,
            extra_count: 0,
        },
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
pub fn filter_strands(sketch: &[KmerCount], ratio_cutoff: f64) -> Vec<KmerCount> {
    let mut filtered = Vec::new();
    for kmer in sketch {
        // "special-case" anything with fewer than 16 kmers -> these are too stochastic to accurately
        // determine if something is an adapter or not. The odds of randomly picking less than 10%
        // (0 or 1 reversed kmers) in 16 should be ~ 17 / 2 ** 16 or 1/4000 so we're avoiding
        // removing "good" kmers
        if kmer.count < 16 {
            filtered.push(kmer.clone());
            continue;
        }

        // check the forward/reverse ratio and only add if it's within bounds
        let lowest_strand_count: u64 = cmp::min(kmer.extra_count, kmer.count - kmer.extra_count);
        if (lowest_strand_count as f64 / kmer.count as f64) >= ratio_cutoff {
            filtered.push(kmer.clone());
        }
    }
    filtered
}

#[test]
fn test_filter_strands() {
    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 10,
            extra_count: 1,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 10,
            extra_count: 2,
        },
        KmerCount {
            hash: 3,
            kmer: vec![],
            count: 10,
            extra_count: 8,
        },
        KmerCount {
            hash: 4,
            kmer: vec![],
            count: 10,
            extra_count: 9,
        },
    ];
    let filtered = filter_strands(&sketch, 0.15);
    assert_eq!(filtered.len(), 4);
    assert_eq!(filtered[0].hash, 1);
    assert_eq!(filtered[3].hash, 4);

    let sketch = vec![
        KmerCount {
            hash: 1,
            kmer: vec![],
            count: 16,
            extra_count: 1,
        },
        KmerCount {
            hash: 2,
            kmer: vec![],
            count: 16,
            extra_count: 2,
        },
        KmerCount {
            hash: 3,
            kmer: vec![],
            count: 16,
            extra_count: 8,
        },
        KmerCount {
            hash: 4,
            kmer: vec![],
            count: 16,
            extra_count: 9,
        },
    ];
    let filtered = filter_strands(&sketch, 0.15);
    assert_eq!(filtered.len(), 2);
    assert_eq!(filtered[0].hash, 3);
    assert_eq!(filtered[1].hash, 4);
}
