# finch #

# Overview #

Finch is an implementation of min-wise independent permutation locality sensitive hashing ("MinHashing") for application in genomic data ("Mashing").
This is a rewrite of One Codex's existing internal clustering/sequence search tool in Rust with a new command-line interface.
We have three primary goals with this rewrite: improve performance, add support for tracking kmer/hash abundances, and improve error filtering.

In addition to some of the memory safety guarentees that Rust enables, we also see considerable speed gains over existing implementations.
Hash generation should hopefully be the rate-limiting step in MinHashing and profiling indicates Finch is spending about a third of its time in `murmurhash3_x64_128` so we should be within an order of magnitude of this theoretical limit.
Some rough benchmarks run on an Early 2015 Macbook Pro: a 4.8 Gb FASTQ (SRR5132341.fastq) takes 99s for Finch to sketch versus 238s for Mash (commit 23776dbe368d398639ec40f133edc06329dc3da8) or 518s for sourmash (commit 5da5ee7c72281ff05cb90d6ce3e8bc4d316998c5).

Finch is one of the first MinHash implementations to support tracking the abundances of each unique kmer/hash.
This allows useful quick interpretations of data (e.g. visualizing sequencing depth with `finch hist`) in addition to laying the support for some of our more advanced filtering features.
Because Finch also supports the [new "Mash" JSON compatibility format](https://github.com/marbl/Mash/blob/master/src/mash/schema-1.0.0.json), we're able to save this count information for downstream processing and comparisons.

 > :warning: Note that although Finch supports the common JSON interchange format, there may still be incompatibilites due to incompaible hashing seed values. :warning:

![Accuracy versus sequencing depth for different filtering schemes](paper/depth_distances.png)
Our default filtering implementation analyses the depth distribution of kmers and dynamically tries to find a cutoff that will remove low abundance (presumably error) kmers as long as their abundance is below a kmer length weighted error rate.
The difference between the "adaptive" approach and a blanket "unique" filter at high sequencing depths of an isolate can be seen in the above figure where many error kmers appear twice or more at moderately high depths (>20) and degrade comparisons against known references.
This approach also balances requirements for metagenomic sequencing by ensuring that only very low abundance organisms whose coverage is similar to the error rate are filtered out (which we find acceptable for sketching comparisons, but isn't an viable strategy for e.g. assembly error trimming).
We also apply a "strand filter" to remove kmers that are only (or close to only) found in one orientation; these are generally sequencing adapters and aren't representative of the underlying samples.

Rather than having to calibrate the size of a Bloom filter or perform exact kmer counting before sketching, we generate a sketch of much larger size than needed (200x) and then filter this "oversketch" into the sketch we need using the histogram of counts contained in the "oversketch".
For the FASTQ above (SRR5132341.fastq), this results in a 60.0 Mb maximum memory footprint for Finch with default filtering and a size 2000 result sketch.
In comparison, Mash has a 1.2 Mb memory footprint and sourmash has a 13.9 Mb memory footprint (both run without filtering options).
While the oversketching approach uses more memory, it's still well within the capabilities of an average laptop or server.


TODO: submit to JOSS (need my orchid ID)
   TODO: run script to generate JOSS metadata - https://gist.github.com/arfon/478b2ed49e11f984d6fb

# Usage #

Finch supports four primary operations, but many of these operations take similar parameters.

## Shared Parameters ##

Finch can take many parameters to control how sketching is performed (and these options are also available for the `dist`, `hist`, and `info` commands).
 - `-n <N>` / `--n-hashes <N>` controls the overall size of the sketch (higher values give better resolution in comparisons).
 - `-k <K>` / `--kmer-length <K>` sets the size of the kmers to be hashed (higher values make comparisons much more taxonomically specific).
 - `--seed <S>` sets the seed for hashing. This should only be changed if directly exporting sketches for comparison with other versions of the Mash algorithm that use a non-zero default seed.

After sketching, filtering is performed and can be controlled through several options:
 - `-f` / `--filter` / `--no-filter` determines whether filtering is applied or not (if not specified, filtering is performed for FASTQ files and not performed for FASTA files by default)
 - `--min-abun-filter <MIN>` / `--max-abun-filter <MAX>` sets absolute minimum and maximum abundances (inclusive) that all kmers must be present at. The minimum filter will override any adaptive error filter guessing (below).
 - `--err-filter <ERR_VALUE>` is the default adaptive filtering scheme. Conceptually `ERR_VALUE` should be approximately the error rate of the sequencer or higher. Note that high error rate sequencing data (i.e. long read sequencing) does poorly with MinHashing in general.
 - `--strand-filter <V>` sets the strand filter cutoff.

If there aren't enough kmers left from the oversketch to satisfy the sketch size, sketching will fail.
There are two options that may help:
 - `--oversketch <N>` can be used to increase the size of the oversketch (normally 200x) and increase the likelihood that the filtered version will be big enough.
 - `--no-strict` will allow finch to proceed with a sketch smaller than the specified size. Use with caution.


## `finch sketch` ##

`finch sketch` will read through a FASTA or FASTQ file and generate a "sketch" that can be used for further .

If the file being read is a sketch and the filtering flag is set (`-f`/`--filter`), the abundance filters will be re-applied to the sketch to allow post-sketch filtering.
The strand filters can only be applied to raw FAST(A/Q) files though, as the strand-level data is lost when sketches are saved.


By default `sketch` generates `.sk` files next to the FAST(A/Q) files it sketches.
This can be overridden by passing either `-O` (capital O) to write to standard out or `-o <file>` (lowercase o) to write to a file.
To read from standard input, use a filename of `-`; this allows streaming of files into `finch`, e.g. `cat testfile.fq | finch sketch -o testfile.sk`.

Sketches should be compatible with the original Mash if you edit their `src/mash/hash.h` and set the hash value to `0` or if you manually override finch's seed value by setting `--seed 42`.

## `finch dist` ##

`finch dist` will calculate Jaccard distances between different sketches.
If the files listed are FAST(A/Q) instead of sketches, Finch will automatically sketch them into memory using the same command-line parameters as in `sketch` or, for files after the first, the parameters used to sketch the first file.

Distances and containments will only be computed from one or more queries to a collection of references; by default the first sketch in the list will be used as the query and all other sketches will be used as references.
This behavior can be manually overriden and other sketches can be used as references by passing `--queries <sketch_1>,<sketch_2`.
Additionally, passing the `--pairwise` option will calculate the distances between all sketches and each other.

Due to different counting algoritms and stopping criteria, distances may be slightly different from the calculation in the original Mash program.
If you'd like comparable (i.e. identical) distances, the `--mash` flag will use their original counting and distance algorithm.

## `finch hist` ##

`finch hist` will output a histogram in JSON format for each sketch provided.
The histogram is a list of the number of hashes at each depth, e.g. `{"sketch_name": [1, 0, 1]}` for a sketch with two hashes at depths 1 and 3 respectively.

## `finch info` ##

`finch info` will output a formatted list of % GC, coverage, etc for each sketch provided.
The values returned from this are approximate and the algoritms used to calculate are still very rough and liable to change.


# References #

There are several other implementations of the Mash algorithm which should be compatible/comparable with this one, notably:
 - [Mash](https://github.com/marbl/Mash) - First implementation and theoretical paper
 - [SourMash](https://github.com/dib-lab/sourmash) - Newer implementation in Python; provides a number of experimental features
