#!/usr/bin/env python

from annotator import Annotator
from annotator import annotate
import gffutils

CHROMS = ['III']
DB_FILE = 'test/data/c_elegans.PRJNA13758.WS257.canonical_geneset.chrIII.25000.gtf.db'
SPECIES = 'ce11'
ANNOTATOR = Annotator.Annotator(db_file=DB_FILE, chroms=CHROMS, species=SPECIES)

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
    with open('test/data/c_elegans.PRJNA13758.WS257.canonical_geneset.H10E21.3a.3UTR.gtf', 'r') as f:
        for line in f:
            utr3_features.append(gffutils.feature.feature_from_line(line))
    return utr3_features

def get_5utr_features_1():
    utr5_features = []
    with open('test/data/c_elegans.PRJNA13758.WS257.canonical_geneset.H10E21.3a.5UTR.gtf', 'r') as f:
        for line in f:
            utr5_features.append(gffutils.feature.feature_from_line(line))
    return utr5_features

def test_get_all_cds_dict_1(a=ANNOTATOR):
    print(
        "Ensures that we're properly setting CDS coordinates (-). "
        "Transcript: H10E21.3a"
        "Start: 12192, End: 14753"
    )
    cds = a.cds_dict['H10E21.3a']
    assert cds['low'] == 12192
    assert cds['hi'] == 14753

def test_annotate_prioritize_cds_1(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (-). "
          "Interval: 'III', '13465', '13466', 'H10E21.3', '0', '-' "
          "Should return CDS"
          )
    chrom = 'III'
    start = 13465
    stop = 13466
    name = 'H10E21.3'
    score = 0
    strand = '-'

    region_priority = [
        ['protein_coding', 'CDS'], ['protein_coding', 'exon'],
    ]

    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    print(annotation)
    assert rname == 'H10E21.3b'
    assert priority == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_cds_2(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (-). "
          "Interval: 'III', '13465', '13466', 'H10E21.3', '0', '-' "
          "Should return CDS even if it's the 2nd priority."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '13465', '13466', 'H10E21.3', '0', '-']
    # )

    chrom = 'III'
    start = 13465
    stop = 13466
    name = 'H10E21.3'
    score = 0
    strand = '-'
    region_priority = [
        ['protein_coding', '3UTR'], ['protein_coding', 'CDS'],
    ]

    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    print(annotation)
    assert rname == 'H10E21.3a' or rname == 'H10E21.3b'
    assert priority == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_1(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'III', '16600', '16601', 'H10E21.1b', '0', '+' "
          "Should return CDS."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '16600', '16601', 'H10E21.1b', '0', '+']
    # )
    region_priority = [
        ['protein_coding', 'CDS'],
        ['protein_coding', 'three_prime_utr'],
        ['protein_coding', 'five_prime_utr'],
    ]
    chrom = 'III'
    start = 16600
    stop = 16601
    name = 'H10E21.1b'
    score = 0
    strand = '+'
    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    print(annotation)
    assert rname == 'H10E21.1b'
    assert priority == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_2(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'III', '16600', '16601', 'H10E21.1b', '0', '+' "
          "Should return 5'UTR."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '16600', '16601', 'H10E21.1a', '0', '+']
    # )
    region_priority = [
        ['protein_coding', 'five_prime_utr'],
        ['protein_coding', 'CDS'],
        ['protein_coding', 'three_prime_utr'],
    ]
    chrom = 'III'
    start = 16600
    stop = 16601
    name = 'H10E21.1b'
    score = 0
    strand = '+'
    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    print(annotation)
    assert rname == 'H10E21.1a'
    assert priority == 'five_prime_utr'
    assert type == 'protein_coding'

def test_annotate_prioritize_3(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'III', '16600', '16601', 'H10E21.1b', '0', '+' "
          "Should return 5'UTR."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '16600', '16601', 'H10E21.1a', '0', '+']
    # )
    region_priority = [
        ['protein_coding', 'three_prime_utr'],
        ['protein_coding', 'five_prime_utr'],
        ['protein_coding', 'CDS'],
    ]
    chrom = 'III'
    start = 16600
    stop = 16601
    name = 'H10E21.1b'
    score = 0
    strand = '+'
    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    print(annotation)
    assert rname == 'H10E21.1a' or rname == 'H10E21.1b'  # don't know which transcript is returned, should clear that up
    assert priority == 'five_prime_utr'
    assert type == 'protein_coding'

def test_annotate_prioritize_4(a=ANNOTATOR):
    print("Tests annotation priority for overlapping transcripts (-/+). "
          "Interval: 'III', '16600', '16601', 'H10E21.1b', '0', '+' "
          "Should return 5'UTR."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '16600', '16601', 'H10E21.1a', '0', '+']
    # )
    chrom = 'III'
    start = 16600
    stop = 16601
    name = 'H10E21.1b'
    score = 0
    strand = '+'
    region_priority = [
        ['protein_coding', 'three_prime_utr'],
        ['protein_coding', 'five_prime_utr'],
        ['protein_coding', 'CDS'],
    ]

    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    gene, rname, priority, type, annotation = a.annotate(
        chrom, start, stop, name, score, strand,
        stranded, region_priority, region_priority
    )
    print(annotation)
    assert rname == 'H10E21.1a' or rname == 'H10E21.1b'  # don't know which transcript is returned, should clear that up
    assert priority == 'five_prime_utr'
    assert type == 'protein_coding'