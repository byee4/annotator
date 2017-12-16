# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)

## [Unreleased]

## [0.0.6] - 2017-12-15
### Added
- Functionality to determine prox vs distal introns (500bp threshold)
- Classify_transcript_type() allow for finer control of annotating noncoding regions.
- Unit tests for proxdist functions.

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

