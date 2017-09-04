#!/usr/bin/env python

from annotator import annotate_bed as a
from annotator import annotate
import gffutils

CHROMS = ['chr19']
DB_FILE = 'test/data/gencode.v19.annotation.chr19.10000000.gtf.db'
SPECIES = 'hg19'

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

def test_get_all_cds_dict_1():
    print(
        "Ensures that we're properly setting CDS coordinates. "
        "Transcript: ENST00000159111.4"
        "Start: 5032902, End: 5151519"
    )
    cds = cds_dict['ENST00000159111.4']
    assert cds['low'] == 5032902
    assert cds['hi'] == 5151519

def test_classify_utr_1():
    print(
        "Tests whether or not we're correctly classifying "
        "three prime utrs given our CDS dictionary and strand (+) info. "
        "Coords: gencode.v19.annotation.ENST00000159111.4.3UTR.gtf"
    )
    utr_features = get_3utr_features_1()
    for feature in utr_features:
        assert a.classify_utr(feature, cds_dict) == '3utr'

def test_classify_utr_2():
    print(
        "Tests whether or not we're correctly classifying "
        "five prime utrs given our CDS dictionary and strand (+) info. "
        "Coords: gencode.v19.annotation.ENST00000159111.4.5UTR.gtf"
    )
    utr_features = get_5utr_features_1()
    for feature in utr_features:
        assert a.classify_utr(feature, cds_dict) == '5utr'


def test_annotate_prioritize_cds_1():
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'chr19', '5071035', '5071036', 'KDM4B', '0', '+' "
          "Should return CDS"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '5071035', '5071036', 'KDM4B', '0', '+']
    # )
    qchrom = 'chr19'
    qstart = 5071035
    qstop = 5071036
    qname = 'KDM4B'
    qscore = 0
    qstrand = '+'
    region_priority = [
        ['protein_coding', 'CDS'], ['protein_coding', 'exon'],
    ]
    stranded = True
    transcript_id_key = 'transcript_id'
    type_key = 'transcript_type'

    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    # assert get_transcript_id(priority) == 'ENST00000159111.4'
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert region == 'CDS'
    assert type == 'protein_coding'

def test_annotate_prioritize_noncoding_exon_1():
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'chr19:5071035-5071036', 'KDM4B', '0', '+' "
          "Should return a retained intron (noncoding exon)"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '5071035', '5071036', 'KDM4B', '0', '+']
    # )
    qchrom = 'chr19'
    qstart = 5071035
    qstop = 5071036
    qname = 'KDM4B'
    qscore = 0
    qstrand = '+'
    region_priority = [
        ['non_coding', 'exon'], ['protein_coding', 'CDS']
    ]

    stranded = True,
    transcript_id_key = 'transcript_id'
    type_key = 'transcript_type'

    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    # assert get_transcript_id(priority) == 'ENST00000592175.1'
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert region == 'noncoding_exon'

def test_annotate_prioritize_noncoding_exon_2():
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
    qchrom = 'chr19'
    qstart = 5071035
    qstop = 5071036
    qname = 'KDM4B'
    qscore = 0
    qstrand = '+'
    stranded = True,
    transcript_id_key = 'transcript_id'
    type_key = 'transcript_type'

    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert region == 'noncoding_exon'

def test_annotate_cds_2():
    print("Tests annotation priority for overlapping transcripts (-). "
          "Interval: 'chr19:4852152-4852191', 'PLIN3', '0', '-' "
          "Should return CDS"
          )
    # interval = pybedtools.create_interval_from_list(
    #     ['chr19', '4852152', '4852191', 'PLIN3', '0', '-']
    # )
    qchrom = 'chr19'
    qstart = 4852152
    qstop = 4852191
    qname = 'PLIN3'
    qscore = 0
    qstrand = '-'
    region_priority = [
        ['protein_coding', 'CDS'], ['non_coding', 'exon'],
    ]
    stranded = True,
    transcript_id_key = 'transcript_id'
    type_key = 'transcript_type'

    chrom, start, end, name, score, strand, \
    gene, rname, region, type, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    # assert get_transcript_id(priority) == 'ENST00000589163.1'
    assert region == 'CDS'
    assert rname == 'PLIN3'
    assert gene == 'ENSG00000105355.4'