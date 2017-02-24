
# finch #

## Overview ##

This is a rewrite of One Codex's internal clustering/sequence search tool in Rust with  new command-line interface.
It also supports the new "Mash" JSON compatibility format (see https://github.com/marbl/Mash/blob/master/src/mash/schema-1.0.0.json) and kmer-abundance counting.

One of our goals is speed. Hash generation should be one of the rate-limiting steps in minhashing and profiling indicates finch is spending about a third of its time in `murmurhash3_x64_128` so we should be within an order of magnitude of theoretical limits. Some rough benchmarks run on a Macbook:

7.5s for finch-rs (commit b23c11fb87992341eaab5d88b53505937f61ba6b)
15.1s for Mash (commit 23776dbe368d398639ec40f133edc06329dc3da8)
~15s for sourmash (-k 21 -n 1000, but only one signature generated?; commit 9b9f7f467b9e99a4bba2e73869058da6021a0835)


## Usage ##

 - `sketch`

 - `dist`

## References ##

 - MASH (https://github.com/marbl/Mash)
 - SourMash (https://github.com/dib-lab/sourmash)
