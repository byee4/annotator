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
python setup.py install
```

### Download Database File (sqlite db files created using [gffutils](https://pythonhosted.org/gffutils/))

[Human hg19 Gencode v19](https://s3-us-west-1.amazonaws.com/genome-references/gencode.v19.annotation.gtf.db)

[Mouse mm10 Gencode vM10](https://s3-us-west-1.amazonaws.com/genome-references/gencode.vM10.annotation.db)

[C. elegans WS257 extended*](https://s3-us-west-1.amazonaws.com/genome-references/c_elegans.PRJNA13758.WS257.canonical_geneset.extend_utr.gtf.db)

* WS257 canonical geneset with extra annotations for annotating upstream/downstream transcripts

In theory any other db file (built from a GTF file) should work... but use at your own risk!!!

# Annotator Example Usage:

```
annotator \
--input BED6_FILE \
--output OUTPUT_FILE \
--gtfdb gencode.v19.annotation.gtf.db \
--species hg19
```

### Output file:

Outputs each original BED interval:

- Chromosome
- Start
- Stop
- Name
- Score
- Strand

Plus annotation stuff (tabbed):

- Assigned Gene ID
- Assigned Gene Name
- Genic Region
- Genic Region Type
- All overlapping annotations


All overlapping annotations are ```|``` seperated, and should follow this format:

```transcript_id:region_start:region_stop:strand:region:gene_id:gene_name:transcript_type:overlap```


### Default Priority (hg19 gencode v19):

- protein_coding CDS
- protein_coding start_codon
- protein_coding stop_codon
- protein_coding 5utr
- protein_coding 3utr
- protein_coding proxintron500
- protein_coding distintron500
- protein_coding Selenocysteine
- non_coding exon
- non_coding proxintron500
- non_coding distintron500
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
randomly appended to the end after all explicitly prioritized features.

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

# create_region_bedfiles
Given a gtf db file, create one or more bedfiles describing merged regions of interest.

###Usage:
```
create_region_bedfiles \
--db_file inputs/gencode.vM10.annotation.db \  # database file downloaded from above
--species mm10 \  # sets the gff/gtf nomenclature (essentially whether it's gencode or wormbase gtf format)
--cds_out outputs/mm10_vM10_cds.bed \  # output cds region
--proxintron_out outputs/mm10_vM10_proxintrons.bed \  # output proximal intron regions (500nt from exons)
--distintron_out outputs/mm10_vM10_distintrons.bed \  # output distal intron regions (> 500nt from exons)
--utr5_out outputs/mm10_vM10_five_prime_utrs.bed \  # output 5'UTR regions
--utr3_out outputs/mm10_vM10_three_prime_utrs.bed  # output 3' UTR regions
```

### Methods:
- To find CDS: parse the gtf file for all 'CDS' featuretypes, then for each gene,
merge overlapping regions.
- To find prox/distal introns: parse the gtf file for all 'exon' featuretypes, then for each
transcript, infer all introns. For each set of introns, classify its distance from an exon as
being either 'proximal' or 'distal' by 500nt (hardcoded for now), and group each region by its
gene id. Then merge overlapping regions on a per-gene basis.
- To find 5/3' UTRs: parse the gtf file for all 'UTR' featuretypes, then use
the CDS features to classify each according to whether or not it lies upstream or
downstream of a transcript's CDS. Then merge overlapping regions on a per-gene basis.

### Notes:
- You do NOT need to specify all output files, just the regions you are interested in
- However in order for this to work with clip_analysis_legacy/analyze_motifs, you DO need to name the outputs as:
    - ${SPECIES}_cds.bed
    - ${SPECIES}_distintron500.bed
    - ${SPECIES}_proxintron500.bed
    - ${SPECIES}_three_prime_utrs.bed
    - ${SPECIES}_five_prime_utrs.bed

# miRNA_name2id:
This script takes a file containing a column with miRNA names and appends an appropriate accession ID
given either a custom name -> accession file, or a gffdb file.

### Usage:
```
miRNA_name2id.py \
--input inputs/all_mirnas.csv \  # tab or SEP separated file
--sep , \  # specify whether or not this file is tab, comma, or some other separated
--custom inputs/ensembl2name_GCm38_mart_export.tsv \  # if applicable, use a custom file (see below)
--name_col miRNA \  # identify the column where the name field is held
--output outputs/all_mirnas.with_ids.csv \
--gffdb inputs/mmu.mirBase_v21.GCm38.gff3.db  # the database file created from a gff file downloaded from mirbase.

```
### Other Options:
```--add_mature``` species whether or not you want an additional column (mature) to your table.

### Using custom versus GFF files for translation:
This script can take either (or both) a gffdb file that you can download from mirbase, or a tabbed custom file
from somewhere else (like ensembl biomart). This file is expected to have the format ```PRECURSOR_ACCESSION\tNAME\tMATURE```
where MATURE is an optional column linking the precursor accession to its processed mature transcript.


### Notes:
- If a miRNA name is associated with more than one accession, the script will return all accessions delimited with ```|```
- You can use the ```build_gffutils_db``` script to create a gffdb from a gff file downloaded from mirbase.

# General notes:
- The ```--species``` flag is only important for setting the GTF file nomenclature; different GTF/GFF files have
differently formatted "attributes" terminologies (see the difference between a wormbase.org and a gencode GTF file).
Currently setting this flag to either 'ce10' or 'ce11' will assuming it is formatted the wormbase way. Otherwise
the format will default to the 'gencode' nomenclature. To use another format, you will have to look at annotate_bed.get_keys(),
which tells this package what keys to expect from the 'attributes' section of the GTF file.
-