use minhashes::KmerCount;
use serialization::SketchDistance;


pub fn distance(sketch1: &Vec<KmerCount>, sketch2: &Vec<KmerCount>, sketch1_name: &str, sketch2_name: &str, mash_mode: bool) -> Result<SketchDistance, &'static str> {
    if sketch1[0].kmer.len() != sketch2[0].kmer.len() {
        return Err("Sketches have different sized kmers");
    }
    let distances;
    if mash_mode {
        distances = calc_distance_mash(sketch1, sketch2);
    } else {
        distances = calc_distance(sketch1, sketch2);
    }
    let containment = distances.0;
    let jaccard = distances.1;
    let common = distances.2;
    let total = distances.3;
    let mash_distance: f64 = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / sketch1[0].kmer.len() as f64;
    Ok(SketchDistance {
        containment: containment,
        jaccard: jaccard,
        mashDistance: f64::min(1f64, f64::max(0f64, mash_distance)),
        commonHashes: common,
        totalHashes: total,
        query: sketch1_name.to_string(),
        reference: sketch2_name.to_string(),
    })
}


fn calc_distance_mash(sketch1: &Vec<KmerCount>, sketch2: &Vec<KmerCount>) -> (f64, f64, u64, u64) {
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


fn calc_distance(sketch1: &Vec<KmerCount>, sketch2: &Vec<KmerCount>) -> (f64, f64, u64, u64) {
    let mut j: usize = 0;
    let mut common: u64 = 0;
    let mut total: u64 = 0;
    
    for i in 0..sketch1.len() {
        while (sketch2[j].hash < sketch1[i].hash) && (j < sketch2.len() - 1) {
            j += 1;
        }
        
        if sketch2[j].hash == sketch1[i].hash {
            common += 1;
        }

        total += 1;
    }

    // Numerator is A-intersect-B, |A| is the denominator, we enforce |A| == |B|
    let containment: f64 = common as f64 / total as f64;
    let jaccard: f64 = common as f64 / (common + 2 * (total - common)) as f64;
    (containment, jaccard, common, total)
}
