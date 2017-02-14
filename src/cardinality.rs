use minhashes::{KmerCount};


pub fn cardinality(sketch: &Vec<KmerCount>) -> Result<u64, &'static str> {
    // Other (possibly more accurate) possibilities:
    // "hyper log-log" estimate from lowest value?
    // multiset distribution applied to total count number?
    // "AKMV" approach: http://people.mpi-inf.mpg.de/~rgemulla/publications/beyer07distinct.pdf
    
    // fast and simple k-minimum value estimate
    // https://research.neustar.biz/2012/07/09/sketch-of-the-day-k-minimum-values/
    return (sketch.len() - 1) / (sketch.last() / usize::max_value())
}
