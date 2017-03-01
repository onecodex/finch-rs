
# finch #

## Overview ##

TODO: Advantages
- Speed
- Filtering
- Counts


This is a rewrite of One Codex's internal clustering/sequence search tool in Rust with  new command-line interface.
It also supports the new "Mash" JSON compatibility format (see https://github.com/marbl/Mash/blob/master/src/mash/schema-1.0.0.json) and kmer-abundance counting.

One of our goals is speed. Hash generation should hopefully be the rate-limiting step in minhashing and profiling indicates finch is spending about a third of its time in `murmurhash3_x64_128` (so we should be within an order of magnitude of this theoretical limit). Some rough benchmarks run on a Macbook:

7.5s for finch-rs (commit b23c11fb87992341eaab5d88b53505937f61ba6b)
15.1s for Mash (commit 23776dbe368d398639ec40f133edc06329dc3da8)
~15s for sourmash (-k 21 -n 1000, but only one signature generated?; commit 9b9f7f467b9e99a4bba2e73869058da6021a0835)


## Usage ##

Finch supports four primary operations:

### `finch sketch` ###

TODO: explain output parameters
TODO: explain sketching parameters

By default `sketch` generates `.sk` files next to the FAST(A/Q) files it sketches.

Sketches should be compatible with the original Mash if you edit `src/mash/hash.h` and set the hash value to `0` or manually override finch's seed value by setting `--seed 42`.

### `finch dist` ###

TODO: explain query/pairwise options

If the files listed are FAST(A/Q) instead of sketches, Finch will automatically sketch them into memory using the same command-line parameters as in `sketch` or, for files after the first, the parameters used to sketch the first file.

Due to different counting algoritms and stopping criteria, distances may be slightly different from the original Mash's.
If you'd like identical distances, the `--mash-mode` flag will use the original counting algorithm.

### `finch hist` ###

TODO: explain

### `finch info` ###

TODO: explain

## References ##

There are several other implementations of the Mash algorithm which should be compatible with this one, notably:
 - [Mash](https://github.com/marbl/Mash) - First implementation and theoretical paper
 - [SourMash](https://github.com/dib-lab/sourmash) - Newer implementation in Python; provides a number of experimental features
