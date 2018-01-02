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

MinHash [@Broder1997] is a document similarity estimation technique that has been applied to problems in genomics including sequence search, phylogenetic reconstruction [@Ondov2016; @Brown2016], and evaluating outbreaks of hospital acquired infections (HAIs) [@Sim2017].
We implement two additions to existing MinHash schemes: calculating abundances (i.e., minmer counts) during the generation of the MinHash sketches and adaptively correcting for biases introduced due to variable sequencing depth. 
This count information greatly improves the utility of MinHashing when working directly from raw read data (i.e., FASTQ files) and allows more robust estimation of distances between both isolates and complex metagenomic samples. 

`finch-rs` is a Rust library and a corresponding command line tool, `finch`, for creating and manipulating MinHash sketches with abundance information.

# References
