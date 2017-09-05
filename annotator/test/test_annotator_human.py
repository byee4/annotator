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

def test_find_introns_1():
    print(
        "Tests the basic functionality of the intron inferring "
        "method. Assumes 'start' and 'end' refer to the lowest "
        "and highest coordinates regardless of strand. This"
        "test tests the lowest negatively stranded transcript."
    )
    transcript_id = 'ENST00000394173.4'
    transcript = transcripts_dict[transcript_id]
    exons = exons_dict[transcript_id]
    introns = a.find_introns(transcript, exons)
    assert introns[0]['start'] == 7808127
    assert introns[0]['end'] == 7808992

def test_find_introns_2():
    print(
        "Tests the basic functionality of the intron inferring "
        "method. Assumes 'start' and 'end' refer to the lowest "
        "and highest coordinates regardless of strand. This"
        "test tests the lowest positively stranded transcript."
    )
    transcript_id = 'ENST00000587541.1'
    transcript = transcripts_dict[transcript_id]
    exons = exons_dict[transcript_id]
    introns = a.find_introns(transcript, exons)
    assert introns[0]['start'] == 490040
    assert introns[0]['end'] == 501668

def test_find_introns_3():
    print(
        "Tests the basic functionality of the intron inferring "
        "method. Assumes 'start' and 'end' refer to the lowest "
        "and highest coordinates regardless of strand. This"
        "tests that a single exon transcript has no intron"
    )
    transcript_id = 'ENST00000589943.1'
    transcript = transcripts_dict[transcript_id]
    exons = exons_dict[transcript_id]
    introns = a.find_introns(transcript, exons)
    assert len(introns) == 0

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
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    # assert get_transcript_id(priority) == 'ENST00000159111.4'
    assert rname == 'KDM4B'
    assert gene == 'ENSG00000127663.10'
    assert region == 'CDS'

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
    gene, rname, region, annotation = a.annotate(
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
    gene, rname, region, annotation = a.annotate(
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
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    # assert get_transcript_id(priority) == 'ENST00000589163.1'
    assert region == 'CDS'
    assert rname == 'PLIN3'
    assert gene == 'ENSG00000105355.4'


def test_intergenic_1():
    print("Tests a region that is intergenic (+)")

    qchrom = 'chr19'
    qstart = 10050000
    qstop = 10006000
    qname = 'intergenic'
    qscore = 0
    qstrand = '+'
    region_priority = [
        ['protein_coding', '3utr'],
        ['protein_coding', '5utr'],
        ['protein_coding', 'CDS'],
    ]

    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    print(annotation)
    assert rname == 'intergenic'  # don't know which transcript is returned, should clear that up
    assert region == 'intergenic'

def test_utr_classification_1():
    print("This tests the same region (chr11:70266235-70266302:+ "
          "as test_utr_classification_2, but re-orders the priority "
          "such that 5utr is higher than 5utr. Supposed to return 5utr")
    chroms = ['chr11']
    db = 'test/data/gencode.v19.annotation.chr11.70M-71M.gtf.db'
    species = 'hg19'

    geneid_to_name_dict, exons_dict, \
    transcripts_dict, cds_dict, features_dict, \
    cds_key, utr3_key, utr5_key, utr_key, \
    gene_name_key, transcript_id_key, type_key = a.create_definitions(
        db, chroms, species
    )
    qchrom = 'chr11'
    qstart = 70266235
    qstop = 70266302
    qname = 'CTTN'
    qscore = 0
    qstrand = '+'
    region_priority = [
        ['protein_coding', '5utr'],
        ['protein_coding', '3utr'],
        ['protein_coding', 'CDS'],
    ]

    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    assert region == '5utr'

def test_utr_classification_2():
    print("This tests the same region (chr11:70266235-70266302:+ "
          "as test_utr_classification_1, but re-orders the priority "
          "such that 3utr is higher than 5utr. Supposed to return 3utr")
    chroms = ['chr11']
    db = 'test/data/gencode.v19.annotation.chr11.70M-71M.gtf.db'
    species = 'hg19'

    geneid_to_name_dict, exons_dict, \
    transcripts_dict, cds_dict, features_dict, \
    cds_key, utr3_key, utr5_key, utr_key, \
    gene_name_key, transcript_id_key, type_key = a.create_definitions(
        db, chroms, species
    )
    qchrom = 'chr11'
    qstart = 70266235
    qstop = 70266302
    qname = 'CTTN'
    qscore = 0
    qstrand = '+'
    region_priority = [
        ['protein_coding', '3utr'],
        ['protein_coding', '5utr'],
        ['protein_coding', 'CDS'],
    ]

    stranded = True,
    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        features_dict, cds_dict, transcript_id_key, type_key
    )
    assert region == '3utr'