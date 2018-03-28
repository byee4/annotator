# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)

## [Unreleased]

## [0.0.10] - 2018-03-25
### Added
- create_AS_STRUCTURE creates AS_STRUCTURE files to be used with clipper
- added annotation_functions.py which will eventually contain shared funcs among scripts

## [0.0.9] - 2018-03-20
### Added
- get_region_lengths now reports total lengths in addition to average length of each genomic region

## [0.0.8] - 2018-02-13
### Added
- Transcript-level region functionality to exons

## [0.0.7] - 2018-02-01
### Added
- Transcript-level region functionality to proximal and distal introns
- Transcript-level region functionality to CDS
- Transcript-level region functionality to 3' and 5' UTR regions

## [0.0.6] - 2017-12-15
### Added
- Functionality to determine prox vs distal introns (500bp threshold)
- Classify_transcript_type() allow for finer control of annotating noncoding regions.
- Unit tests for proxdist functions.
- build_gffutils_db
- create_region_bedfiles
- gene_name2id
- miRNA_name2id

### Changed
- datasets/*priority.txt to reflect prox and distal intron priorities

### Deprecated
- is_protein_coding() in favor of classify_transcript_type()

## 0.0.5 - 2017-10-18
### Added
- First sharable commit to github
- Lab slides that show usage
- README now contains examples and default priorities/params

[Unreleased]: https://github.com/byee4/annotator...HEAD

