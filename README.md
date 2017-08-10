# annotator

### Requirements: See the environments.yaml file, or:

```
conda create -n annotator \
python=2.7 \
gffutils=0.8.7.1 \
bedtools=2.26 \
pybedtools=0.7.10 \
tqdm=4.14

source activate annotator

conda install --override-channels \
-c conda-forge bzip2 # fixes some weird c-library issue
```

### Installation:
```
git clone https://github.com/byee4/annotator/
cd annotator
python setup.py build
python setup.py install
```

### Download Database File (sqlite db files created using [gffutils](https://pythonhosted.org/gffutils/))

[Human hg19 Gencode v19](https://s3-us-west-1.amazonaws.com/genome-references/gencode.v19.annotation.gtf.db)

[Mouse mm10 Gencode vM10](https://s3-us-west-1.amazonaws.com/genome-references/gencode.vM10.annotation.db)

[C. elegans WS257 extended*](https://s3-us-west-1.amazonaws.com/genome-references/c_elegans.PRJNA13758.WS257.canonical_geneset.extend_utr.gtf.db)

* WS257 canonical geneset with extra annotations for annotating upstream/downstream transcripts

In theory any other db file (built from a GTF file) should work... but use at your own risk!!!

### Example Usage:

```
annotate-bed \
--input BED6_FILE \
--output OUTPUT_FILE \
--gtfdb gencode.v19.annotation.gtf.db \
```

### Output file:

Outputs each original BED interval, plus annotation stuff:

- Chromosome
- Start
- Stop
- Name
- Score
- Strand
- Assigned Gene ID
- Assigned Gene Name
- Genic Region
- Genic Region Type
- All overlapping annotations

outputs bedfile + first priority annotation + all overlapping annotations
in this format:

```transcript_id:region_start:region_stop:strand:region:gene_id:gene_name:transcript_type:overlap```


### Default Priority (hg19 gencode v19):

- protein_coding CDS
- protein_coding start_codon
- protein_coding stop_codon
- protein_coding THREE_AND_FIVE_PRIME_UTR
- protein_coding 5UTR
- protein_coding 3UTR
- protein_coding intron
- protein_coding Selenocysteine
- non_coding exon
- non_coding intron
- non_coding transcript
- non_coding gene
- non_coding Selenocysteine

### Methods:

This script works first by prioritizing overlapping transcript regions
for each gene, then by prioritizing genic regions among all overlapping
genes that overlap your region of interest. As in, if a feature overlaps
both a 3'UTR of Transcript A and a CDS of Transcript B
(both belonging to Gene X), we decide to report Transcript B:CDS.
If multiple genes overlap your feature (Gene X, Gene Y), we prioritize all
transcripts within each gene first such that each gene contains a 'chosen'
transcript, then prioritize each gene in the same way.

### Other Options:

```--transcript-priority-file``` determines the priority when ordering
transcripts within each gene. It's a comma delimited file
(See the ```data/priority.txt``` for an example) containing both
the feature type, transcript type prioritized by line order. Features
that are discovered to be overlapped but are not in this list will be
randomly appended to the end as last priority.

```--gene-priority-file``` determines the priority when choosing which gene
to report. The format is identical to transcript-priority-file.

```--species``` specifies whether or not the database is formatted
to gencodegenes specifications (ie. mm10/hg19) or to wormbase
spec (ce11). Default is hg19/gencodegenes GTF convention.

```--unstranded``` will allow for unstranded features to be selected,
however if strand is specified in the BED file, we'll try to look for
correctly stranded features first. If strand is not specified, we'll
prioritize positive stranded features first.

```--limit-chroms-to``` will limit the dictionary build to only
include these chromosomes for faster processing and less memory
footprint. Leave blank to hash all chromosomes in the db file

Let me know if you have issues/questions: bay001@ucsd.edu
