@0x8c29b313fdc31ca5;

# right now all SketchMethods assume canonicalized, case-insensitive nucleotide records
enum SketchMethod {
  murmurHash3 @0;
  murmurHash3Scaled @1;
  none @2;
}

struct FilterParams {
  filtered @0 :Bool;
  lowAbunFilter @1 :UInt32;
  highAbunFilter @2 :UInt32;
  errFilter @3 :Float64;
  strandFilter @4 :Float64;
}

struct SketchParams {
  sketchMethod @0 :SketchMethod;
  kmerLength @1 :UInt8;
  # common hash-based sketch parameters
  kmersToSketch @2 :UInt64;
  hashSeed @3 :UInt64;
  # parameters for Finch's "mash" sketches
  finalSize @4 :UInt64;
  noStrict @5 :Bool;
  # parameter for scaled sketching
  scale @6 :Float64;
}

# a kmer; the basic unit of the sketch
# note that we don't track k-mer locations because we could potentially have
# to store >1000 positions in here and that's tricky
struct KmerCount {
  hash @0 :UInt64;
  kmer @1 :Data; 
  count @2 :UInt32;
  extraCount @3 :UInt32;
  label @4 :Data;
}

struct Sketch {
  name @0 :Text;
  # useful metadata
  seqLength @1 :UInt64;
  numValidKmers @2 :UInt64;
  comment @3 :Text;
  
  hashes @4 :List(KmerCount);
  filterParams @5 :FilterParams;
  sketchParams @6 :SketchParams;
}

struct Multisketch {
  sketches @0 :List(Sketch);
}
