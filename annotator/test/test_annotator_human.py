#!/usr/bin/env python

from annotator import Annotator
from annotator import annotate
import gffutils

CHROMS = ['chr19']
DB_FILE = 'test/data/gencode.v19.annotation.chr19.10000000.gtf.db'
SPECIES = 'hg19'
ANNOTATOR = Annotator.Annotator(
    db_file=DB_FILE, chroms=CHROMS, species=SPECIES
)

def get_transcript_id(st):
    return st.split(':')[0]

def get_featuretype(st):
    return st.split(':')[4]

def gene_priority_1():
    return annotate.parse('priority.txt')

def transcript_priority_1():
    return annotate.parse('priority.txt')

def get_3utr_features_1():
    utr3_features = []
    with open('test/data/gencode.v19.annotation.ENST00000159111.4.3UTR.gtf', 'r') as f:
        for line in f:
            utr3_features.append(gffutils.feature.feature_from_line(line))
    return utr3_features

def get_5utr_features_1():
    utr5_features = []
    with open('test/data/gencode.v19.annotation.ENST00000159111.4.5UTR.gtf', 'r') as f:
        for line in f:
            utr5_features.append(gffutils.feature.feature_from_line(line))
    return utr5_features

def test_get_all_cds_dict_1(a=ANNOTATOR):
    print(
        "Ensures that we're properly setting CDS coordinates. "
        "Transcript: ENST00000159111.4"
        "Start: 5032902, End: 5151519"
    )
    cds = a.cds_dict['ENST00000159111.4']
    assert cds['low'] == 5032902
    assert cds['hi'] == 5151519

def test_classify_utr_1(a=ANNOTATOR):
    print(
        "Tests whether or not we're correctly classifying "
        "three prime utrs given our CDS dictionary and strand (+) info. "
        "Coords: gencode.v19.annotation.ENST00000159111.4.3UTR.gtf"
    )
    utr_features = get_3utr_features_1()
    for feature in utr_features:
        assert a._classify_utr(feature) == '3UTR'

def test_classify_utr_2(a=ANNOTATOR):
    print(
        "Tests whether or not we're correctly classifying "
        "five prime utrs given our CDS dictionary and strand (+) info. "
        "Coords: gencode.v19.annotation.ENST00000159111.4.5UTR.gtf"
    )
    utr_features = get_5utr_features_1()
    for feature in utr_features:
        assert a._classify_utr(feature) == '5UTR'


def test_annotate_prioritize_cds_1(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'chr19', '5071035', '5071036', 'KDM4B', '0', '+' "
          "Should return CDS"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '5071035', '5071036', 'KDM4B', '0', '+']
    # )
    chrom = 'chr19'
    start = 5071035
    stop = 5071036
    name = 'KDM4B'
    score = 0
    strand = '+'
    region_priority = [
        ['protein_coding', 'CDS'], ['protein_coding', 'exon'],
    ]

    stranded = True,
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    # assert get_transcript_id(priority) == 'ENST00000159111.4'
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert priority == 'CDS'

def test_annotate_prioritize_noncoding_exon_1(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'chr19:5071035-5071036', 'KDM4B', '0', '+' "
          "Should return a retained intron (noncoding exon)"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '5071035', '5071036', 'KDM4B', '0', '+']
    # )
    chrom = 'chr19'
    start = 5071035
    stop = 5071036
    name = 'KDM4B'
    score = 0
    strand = '+'
    region_priority = [
        ['non_coding', 'exon'], ['protein_coding', 'CDS']
    ]

    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    # assert get_transcript_id(priority) == 'ENST00000592175.1'
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert priority == 'noncoding_exon'

def test_annotate_prioritize_noncoding_exon_2(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'chr19:5071035-5071036', 'KDM4B', '0', '+' "
          "Should return a retained intron (noncoding exon)"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '5071035', '5071036', 'KDM4B', '0', '+']
    # )
    region_priority = [
        ['non_coding', 'exon'], ['protein_coding', 'CDS']
    ]
    chrom = 'chr19'
    start = 5071035
    stop = 5071036
    name = 'KDM4B'
    score = 0
    strand = '+'
    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert priority == 'noncoding_exon'

def test_annotate_cds_2(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (-). "
          "Interval: 'chr19:4852152-4852191', 'PLIN3', '0', '-' "
          "Should return CDS"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '4852152', '4852191', 'PLIN3', '0', '-']
    # )
    chrom = 'chr19'
    start = 4852152
    stop = 4852191
    name = 'PLIN3'
    score = 0
    strand = '-'
    region_priority = [
        ['protein_coding', 'CDS'], ['non_coding', 'exon'],
    ]
    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    # assert get_transcript_id(priority) == 'ENST00000589163.1'
    assert priority == 'CDS'
    assert rname == 'PLIN3'
    assert gene == 'ENSG00000105355.4'