---
title: 'Finch'
tags:
  - minhash
  - Rust
authors:
 - name: Roderick Bovee
   orcid: 0000-0002-8819-9549
   affiliation: 1
 - name: Nick Greenfield
   orcid: 0000-0001-8637-406X
   affiliation: 1
affiliations:
 - name: One Codex, San Francisco, California, USA
   index: 1
date: 24 August 2017
bibliography: paper.bib
---
# Summary


MinHash [@Broder1997] is a document similarity estimation technique that has been applied to problems in genomics including sequence search, phylogenetic reconstruction [@Ondov2016; @Brown2016], and evaluating outbreaks of hospital acquired infections (HAIs) [@Sim2017]. We developed the `finch-rs` library (https://github.com/onecodex/finch-rs) and `finch` command line tool for creating, filtering, and manipulating MinHash sketches from genomics data, including both FASTA sequence files and FASTQ raw read data from next-generation sequencing (NGS) instruments. We extend existing MinHash schemes for genomics data with two major additions: (1) calculation of abundances (i.e., minmer counts) during the generation of the MinHash sketches; and (2) adaptive correction of biases introduced due to variable sequencing depths. These features greatly improve the utility of MinHashing when applied directly to raw read data (i.e., FASTQ files) and allows more robust estimation between both isolates and complex metagenomic samples.

Finch and similar genomic MinHashing software works by breaking sequence data up into k-length nucleotide or amino acid subsequences ("k-mers"), computing the hash of each k-mer, and then taking the _n_ lowest hash values. Collectively, these _n_ smallest values ("minmers") comprise a "sketch" of the input sample. By default, previous MinHash implementations for genomics data work by creating sketches from _all_ k-mers from an input genomic dataset (though the original Mash tool does enable filtering out k-mers that appear only once using a Bloom filter [@Ondov2016]). While this works well for high-quality sequences such as genome assemblies (i.e., FASTA files), it quickly becomes problematic when working with raw FASTQ data, where errors from NGS instruments can lead to a far larger number of _unique_ observed k-mers than are truly present biologically. Similarly, this also leads to the inclusion of sequencing errors and k-mers from minor community members when comparing complex, mixed genomic samples (i.e., microbiome samples). In both cases, non-representative k-mers (either direct products of sequencing error or low abundance organisms) come to dominate sketches and confound inter-sample distance estimates.

We address this filtering challenge by "over-sketching" the input genomic data. First, we create a sketch substantially larger than the desired final size (_n_), tracking the abundances of each k-mer in the sketch, and using the abundances in the large sketch to determine a dynamic filtering threshold. We also track how often each k-mer is seen in its forward versus its reverse orientation. These two metrics allow us to both: (1) estimate the empirical sequencing error in the sample and only select k-mers that appear to be biologically present; and (2) remove k-mers that exhibit unbalanced forward and reverse orientation ratios. The former addresses the challenges of comparing data sequenced to varying depths (or samples of varying natural complexity), while the latter can correct for errors that may stem from measurement artifacts such as reads that include adapters and barcode sequences. By removing these error k-mers from a final, reduced sketch (size _n_), `finch` more robustly estimates distances between sets of both isolates and complex metagenomic samples.

Finch is written in the Rust programming language [@Matsakis2014] which reduces programming errors through static type checking, allows greater control over performance, and easily integrates into higher-level languages such as Python or R. We have integrated Finch into our One Codex platform [@Minot2015] and are using it to power clustering and similarity search features.


# References
