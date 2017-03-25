---
title: 'Finch'
tags:
  - minhash
  - Rust
authors:
 - name: Roderick Bovee
   orcid: FIXME
   affiliation: 1
 - name: Nick Greenfield
   orcid: FIXME
   affiliation: 1
affiliations:
 - name: One Codex
   index: 1
date: 23 March 2017
bibliography: paper.bib
---
# Summary

MinHash [@Broder1997] is a document similarity estimation technique that has been applied to the theoretical bioinformatics problems of sequence search and genome comparisons [@Ondov2016; @Brown2016] and used for testing infectious disease transmission hypotheses [@Sim2017].
We implement two additions to existing minhash schemes: depth tracking of the minhashes and filtering based upon these depths.
Depth tracking allows either rough heuristics on either isolate purity or metagenome complexity.

finch contains a Rust library and a command line tool for creating and manipulating MinHash sketches.

# References
