# Finch #

[![CircleCI](https://circleci.com/gh/onecodex/finch-rs.svg?style=shield&circle-token=4fa4296a595f0d0e8dfcea750d76c39576f2df73)](https://circleci.com/gh/onecodex/finch-rs)

Finch is an implementation of min-wise independent permutation locality sensitive hashing ("MinHashing") for application to genomic data ("Mashing").
This repository provides a library and command-line interface that reimplements much of One Codex's [existing internal clustering/sequence search tool](http://blog.onecodex.com/2016/02/11/clustering-and-basespace/) (and adds new features/extensions!) in Rust.

## Design goals
We have 3 primary design goals with Finch:

1. Support for tracking k-mer/minmer _counts_;

2. Improved error filtering for use with raw sequence data (i.e., FASTQs); and

3. Improved performance.


### Support for counts
Finch is one of the first MinHash implementations to support tracking the abundances of each unique _k_-mer/minmer.
This allows useful quick interpretations of data (e.g. visualizing sequencing depth with `finch hist`) in addition to laying the support for more advanced filtering (implemented here) and distance measures (a future area for investigation and study).
Because Finch also supports the [new "Mash" JSON compatibility format](https://github.com/marbl/Mash/blob/master/src/mash/schema-1.0.0.json), we're able to save this count information for downstream processing and comparisons.

 > :warning: &nbsp;Note that although Finch supports the common JSON interchange format, there may still be incompatibilites due to incompaible hashing seed values (Finch uses a seed of `0`, while Mash uses `42`. See [this issue](https://github.com/marbl/Mash/issues/12) for more discussion).

### Filtering
A key design goal for Finch is to stably, accurately sketch variable depth sequencing of bacterial isolates (i.e., to not require a high quality assembly). While Mash includes the ability to specify a fixed absolute cutoff (or use a Bloom filter to filter out unique minmers in a sample), in practice this filtering mechanism performs poorly across varied sequencing depths, genome sizes, and sequencing error profiles. This filtering approach extends several of the ideas in the original Mash paper, as well as [prior work on clustering clinical _C. diff_ isolates](http://www.dmidjournal.com/article/S0732-8893(16)30322-4/abstract).

Finch includes two classes of automated filtering (on by default for FASTQs) to address this challenge:

1. *Adapative, count-based filtering*: Our default implementation "over-sketches" an input file (i.e., includes _n_ fold more minmers than specified by the user) and dynamically determines an appropriate count-based filtering cutoff that removes low abundance _k_-mers (which predominantly consist of random errors). As a fallback, this implementation never excludes more than an upper bound proportion of minmers that may be reasonably assumed to be erroneous (`k` times a specified `--err-filter` error rate, for a default of 21%).

2. *Strandedness filtering*: We also apply a "strandedness filter" to remove _k_-mers that are overwhelmingly found in one orientation; in practice these are often seqeuencing adapters or other artifacts that are not representative of the underlying sample (and may erroneously reduce the similarity between two sketches).

In practice, these optimizations dramatically increase the stability and accuracy of distance estimates between variable depth isolates and finished assemblies. The below graph shows the distance between a high quality bacterial assembly and sketches of FASTQ files representing various sequencing depths across 3 filter strategies: no filtering (`cutoff=0`), "naïve" filtering (`cutoff=1`), and Finch's adapative, variable-cutoff filtering. We use the containment measure, with results closer to `100%` being better:

![Accuracy versus sequencing depth for different filtering schemes](https://github.com/onecodex/finch-rs/blob/master/paper/depth_distance.png?raw=true)

Notably, this strategy performs well across sequencing depths, while fixed cutoffs of 0 and 1 never achieve the near-100% expected containment (and the former fails catastrophically at even modest depths). It is also robust to repeated errors (i.e., error _k_-mers/minmers with a count >= 2), which can begin to pose a problem at sequencing depths as low a 10-20X.

### Speed/performance
In addition to some of the memory safety guarentees that Rust enables, we also see considerable speed gains over existing implementations.
Ideally, hashing should be the rate-limiting step in MinHashing. Profiling indicates `finch` spends about a third of its time in `murmurhash3_x64_128` so we should be within an order of magnitude of this theoretical limit.

|                 | Mash    | Sourmash | Finch (no filtering) | Finch (automatic filtering) |
|-----------------|---------|----------|----------------------|-----------------------------|
| Time (s)        | 238     | 518      | **99**               | 104                         |
| Max Memory (Mb) | **1.2** | 13.9     | 21.8                 | 60.0                        |

> *Note: Benchmarks run on an Early 2015 Macbook Pro. Benchmark is sketching a 4.8Gb FASTQ file ([SRR5132341.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5132341)). Benchmarks used the following commits for Mash ([23776db](https://github.com/marbl/Mash/commit/23776dbe368d398639ec40f133edc06329dc3da8)) and sourmash ([5da5ee7](https://github.com/dib-lab/sourmash/commit/5da5ee7c72281ff05cb90d6ce3e8bc4d316998c5))*

## Usage ##

Finch supports four primary operations, with many of these operations taking similar parameters.

### Shared Parameters ###

Finch can take many parameters to control how sketching is performed (and these options are also available for the `dist`, `hist`, and `info` commands).
 - `-n <N>` / `--n-hashes <N>` controls the overall size of the sketch (higher values give better resolution in comparisons).
 - `-k <K>` / `--kmer-length <K>` sets the size of the kmers to be hashed (higher values make comparisons much more taxonomically specific).
 - `--seed <S>` sets the seed for hashing. This should only be changed if directly exporting sketches for comparison with other versions of the Mash algorithm that use a non-zero default seed.

After sketching, filtering is performed and can be controlled through several options:
 - `-f` / `--filter` / `--no-filter` determines whether filtering is applied or not (if not specified, filtering is performed for FASTQ files and not performed for FASTA files by default)
 - `--min-abun-filter <MIN>` / `--max-abun-filter <MAX>` sets absolute minimum and maximum abundances (inclusive) that all kmers must be present at. The minimum filter will override any adaptive error filter guessing (below).
 - `--err-filter <ERR_VALUE>` is the default adaptive filtering scheme. Conceptually `ERR_VALUE` should be approximately the error rate of the sequencer or higher. Note that high error rate sequencing data (i.e. long read sequencing) does poorly with MinHashing in general.
 - `--strand-filter <V>` sets the strand filter cutoff.

Note that if there aren't enough kmers left from the oversketch to satisfy the sketch size, sketching will fail.
There are two options that may help:
 - `--oversketch <N>` can be used to increase the size of the oversketch (normally 200x) and increase the likelihood that the filtered version will be big enough.
 - `--no-strict` will allow finch to proceed with a sketch smaller than the specified size. Use with caution.


### `finch sketch` ###

`finch sketch` will read through a FASTA or FASTQ file and generate a "sketch" that can be used for further .

If the file being read is a sketch and the filtering flag is set (`-f`/`--filter`), the abundance filters will be re-applied to the sketch to allow post-sketch filtering.
The strand filters can only be applied to raw FAST(A/Q) files though, as the strand-level data is lost when sketches are saved.

By default `sketch` generates `.sk` suffixed files next to the FAST(A/Q) files it sketches.
This can be overridden by passing either `-O` (capital O) to write to standard out or `-o <file>` (lowercase o) to write to a file.
To read from standard input, use a filename of `-`; this allows streaming of files into `finch`, e.g. `cat testfile.fq | finch sketch -o testfile.sk -`.

Sketches should be compatible with the original Mash implementation if you edit their `src/mash/hash.h` and set the hash value to `0` or if you manually override Finch's seed value by setting `--seed 42`.

### `finch dist` ###

`finch dist` will calculate Jaccard distances between different sketches.
If the files listed are FASTA/Q instead of sketches, Finch will automatically sketch them into memory using the same command-line parameters as in `sketch` or, for files after the first, the parameters used to sketch the first file.

Distances and containments will only be computed from one or more queries to a collection of references; by default the first sketch in the list will be used as the query and all other sketches will be used as references.
This behavior can be manually overriden and other sketches can be used as references by passing `--queries <sketch_1>,<sketch_2`.
Additionally, passing the `--pairwise` option will calculate the distances between all sketches and each other.

Due to different counting algoritms and stopping criteria, distances may be slightly different from the calculation in the original Mash program.
If you'd like identical distances, the `--mash` flag will use their original counting and distance algorithm (differences relate to where to stop in the minmer set comparisons).

### `finch hist` ###

`finch hist` will output a histogram in JSON format for each sketch provided.
The histogram is a list of the number of minmers at each depth, e.g. `{"sketch_name": [1, 0, 1]}` for a sketch with two minmers, one with a depth of 1 (first position) and one with a depth of 3 (third position).

### `finch info` ###

`finch info` will output a formatted list of % GC, coverage, etc for each sketch provided.

> :warning: &nbsp;Note that the values returned from this are approximate and the algoritms used to calculate are still rough and liable to change.


## References ##

There are several other implementations of the Mash algorithm which should be compatible/comparable with this one, notably:
 - [Mash](https://github.com/marbl/Mash) - First implementation and theoretical paper
 - [SourMash](https://github.com/dib-lab/sourmash) - Newer implementation in Python; provides a number of experimental features
