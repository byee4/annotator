#!/usr/bin/env python

from annotator import annotate_bed as a
from annotator import annotate
import gffutils

CHROMS = ['III']
DB_FILE = 'test/data/c_elegans.PRJNA13758.WS257.canonical_geneset.chrIII.25000.gtf.db'
SPECIES = 'ce11'

geneid_to_name_dict, exons_dict, \
transcripts_dict, cds_dict, features_dict, \
cds_key, utr3_key, utr5_key, utr_key, \
gene_name_key, transcript_id_key, type_key = a.create_definitions(
    DB_FILE, CHROMS, SPECIES
)

def get_transcript_id(st):
    return st.split(':')[0]

def get_featuretype(st):
    return st.split(':')[4]

def gene_priority_1():
    return annotate.parse_annotation_priority('priority.txt')

def transcript_priority_1():
    return annotate.parse_annotation_priority('priority.txt')

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

def test_get_all_cds_dict_1():
    print(
        "Ensures that we're properly setting CDS coordinates (-). "
        "Transcript: H10E21.3a"
        "Start: 12192, End: 14753"
    )
    cds = cds_dict['H10E21.3a']
    assert cds['low'] == 12192
    assert cds['hi'] == 14753

def test_annotate_prioritize_cds_1():
    print("Tests annotation priority for overlapping transcripts (-). "
          "Interval: 'III', '13465', '13466', 'H10E21.3', '0', '-' "
          "Should return CDS"
          )
    qchrom = 'III'
    qstart = 13465
    qstop = 13466
    qname = 'H10E21.3'
    qscore = 0
    qstrand = '-'

    region_priority = [
        ['protein_coding', 'CDS'], ['protein_coding', 'exon'],
    ]

    stranded = True,
    # transcript_priority = transcript_priority_1()
    # gene_priority = gene_priority_1()
    transcript_id_key = 'transcript_id'
    type_key = 'transcript_biotype'

    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    assert rname == 'H10E21.3b'
    assert region == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_cds_2():
    print("Tests annotation priority for overlapping transcripts (-). "
          "Interval: 'III', '13465', '13466', 'H10E21.3', '0', '-' "
          "Should return CDS even if it's the 2nd priority."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '13465', '13466', 'H10E21.3', '0', '-']
    # )

    qchrom = 'III'
    qstart = 13465
    qstop = 13466
    qname = 'H10E21.3'
    qscore = 0
    qstrand = '-'
    region_priority = [
        ['protein_coding', '3UTR'], ['protein_coding', 'CDS'],
    ]

    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    assert rname == 'H10E21.3a' or rname == 'H10E21.3b'
    assert region == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_1():
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
    qchrom = 'III'
    qstart = 16600
    qstop = 16601
    qname = 'H10E21.1b'
    qscore = 0
    qstrand = '+'
    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    print(annotation)
    assert rname == 'H10E21.1b'
    assert region == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_2():
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
    qchrom = 'III'
    qstart = 16600
    qstop = 16601
    qname = 'H10E21.1b'
    qscore = 0
    qstrand = '+'
    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    print(annotation)
    assert rname == 'H10E21.1a'
    assert region == 'five_prime_utr'
    assert type == 'protein_coding'

def test_annotate_prioritize_3():
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
    qchrom = 'III'
    qstart = 16600
    qstop = 16601
    qname = 'H10E21.1b'
    qscore = 0
    qstrand = '+'
    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    print(annotation)
    assert rname == 'H10E21.1a' or rname == 'H10E21.1b'  # don't know which transcript is returned, should clear that up
    assert region == 'five_prime_utr'
    assert type == 'protein_coding'

def test_annotate_prioritize_4():
    print("Tests annotation priority for overlapping transcripts (-/+). "
          "Interval: 'III', '16600', '16601', 'H10E21.1b', '0', '+' "
          "Should return 5'UTR."
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['III', '16600', '16601', 'H10E21.1a', '0', '+']
    # )
    qchrom = 'III'
    qstart = 16600
    qstop = 16601
    qname = 'H10E21.1b'
    qscore = 0
    qstrand = '+'
    region_priority = [
        ['protein_coding', 'three_prime_utr'],
        ['protein_coding', 'five_prime_utr'],
        ['protein_coding', 'CDS'],
    ]

    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    print(annotation)
    assert rname == 'H10E21.1a' or rname == 'H10E21.1b'  # don't know which transcript is returned, should clear that up
    assert region == 'five_prime_utr'
    assert type == 'protein_coding'