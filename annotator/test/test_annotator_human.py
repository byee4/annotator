#!/usr/bin/env python

from annotator import annotate_bed as a
from annotator import annotation_functions as af
from annotator import annotator
import gffutils
import pybedtools

CHR19_CHROMS = ['chr19']
CHR21_CHROMS = ['chr21']
CHR19_DB_FILE = 'test/data/gencode.v19.annotation.chr19.10000000.gtf.db'
CHR21_DB_FILE = 'test/data/gencode.v19.annotation.chr21.46M-47M.gtf.db'
SPECIES = 'hg19'

chr19_exons_dict, \
chr19_transcripts_dict, chr19_cds_dict, chr19_features_dict, chr19_keys = a.create_definitions(
    CHR19_DB_FILE, CHR19_CHROMS, SPECIES
)

chr21_exons_dict, \
chr21_transcripts_dict, chr21_cds_dict, chr21_features_dict, chr21_keys = a.create_definitions(
    CHR21_DB_FILE, CHR21_CHROMS, SPECIES
)


chr19_transcript_id_key = chr19_keys['transcript_id']
chr19_type_key = chr19_keys['transcript_type']
chr19_utr_key = chr19_keys['utr']

chr21_transcript_id_key = chr21_keys['transcript_id']
chr21_type_key = chr21_keys['transcript_type']
chr21_utr_key = chr21_keys['utr']

def get_transcript_id(st):
    return st.split(':')[0]

def get_featuretype(st):
    return st.split(':')[4]

def gene_priority_1():
    return annotator.parse_annotation_priority('priority.txt')

def transcript_priority_1():
    return annotator.parse_annotation_priority('priority.txt')

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

def get_proxdist_db():
    return 'test/data/proxdist_test.gtf.db'

def test_get_all_cds_dict_1():
    print(
        "Ensures that we're properly setting CDS coordinates. "
        "Transcript: ENST00000159111.4"
        "Start: 5032902, End: 5151519"
    )
    cds = chr19_cds_dict['ENST00000159111.4']
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
    transcript = chr19_transcripts_dict[transcript_id]
    exons = chr19_exons_dict[transcript_id]
    introns = af.find_introns(transcript, exons)
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
    transcript = chr19_transcripts_dict[transcript_id]
    exons = chr19_exons_dict[transcript_id]
    introns = af.find_introns(transcript, exons)
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
    transcript = chr19_transcripts_dict[transcript_id]
    exons = chr19_exons_dict[transcript_id]
    introns = af.find_introns(transcript, exons)
    assert len(introns) == 0

def test_classify_utr_1():
    print(
        "Tests whether or not we're correctly classifying "
        "three prime utrs given our CDS dictionary and strand (+) info. "
        "Coords: gencode.v19.annotation.ENST00000159111.4.3UTR.gtf"
    )
    utr_features = get_3utr_features_1()
    for feature in utr_features:
        assert af.classify_utr(feature, chr19_cds_dict) == '3utr'

def test_classify_utr_2():
    print(
        "Tests whether or not we're correctly classifying "
        "five prime utrs given our CDS dictionary and strand (+) info. "
        "Coords: gencode.v19.annotation.ENST00000159111.4.5UTR.gtf"
    )
    utr_features = get_5utr_features_1()
    for feature in utr_features:
        assert af.classify_utr(feature, chr19_cds_dict) == '5utr'


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

    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        chr19_features_dict, chr19_cds_dict, chr19_keys
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

    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        chr19_features_dict, chr19_cds_dict, chr19_keys
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

    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        chr19_features_dict, chr19_cds_dict, chr19_keys
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

    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, region_priority, region_priority,
        chr19_features_dict, chr19_cds_dict, chr19_keys
    )
    # assert get_transcript_id(priority) == 'ENST00000589163.1'
    assert region == 'CDS'
    assert rname == 'PLIN3'
    assert gene == 'ENSG00000105355.4'


def test_annotate_distintron_1():
    print("Tests annotation priority for overlapping transcripts (+). "
          "Interval: 'chr21:46553734-46553788', 'ADARB1', '0', '+' "
          "Should return distintron"
          )

    qchrom = 'chr21'
    qstart = 46553734
    qstop = 46553788
    qname = 'ADARB1'
    qscore = 0
    qstrand = '+'
    transcript_priority = [
        ['protein_coding', 'distintron500'], ['non_coding', 'proxintron500'],
    ]
    gene_priority = [
        ['non_coding', 'proxintron500'], ['protein_coding', 'distintron500'],
    ]


    stranded = True,

    chrom, start, end, name, score, strand, \
    gene, rname, region, annotation = a.annotate(
        qchrom, qstart, qstop, qname, qscore, qstrand,
        stranded, transcript_priority, gene_priority,
        chr21_features_dict, chr21_cds_dict, chr21_keys
    )
    assert rname == 'ADARB1'
    assert region == 'distintron500'

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
        chr19_features_dict, chr19_cds_dict, chr19_keys
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

    exons_dict, \
    transcripts_dict, cds_dict, features_dict, \
    keys = a.create_definitions(
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
        features_dict, cds_dict, keys
    )
    assert region == '5utr'

def test_utr_classification_2():
    print("This tests the same region (chr11:70266235-70266302:+ "
          "as test_utr_classification_1, but re-orders the priority "
          "such that 3utr is higher than 5utr. Supposed to return 3utr")
    chroms = ['chr11']
    db = 'test/data/gencode.v19.annotation.chr11.70M-71M.gtf.db'
    species = 'hg19'

    exons_dict, \
    transcripts_dict, cds_dict, features_dict, \
    keys = a.create_definitions(
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
        features_dict, cds_dict, keys
    )
    assert region == '3utr'

def test_split_prox_dist_1():
    print("Tests the core functionality of assigning proximal and distal "
          "intron spaces. Should not return any distal introns based on the "
          "specified distance.")
    length = 10
    midpoint = 5
    intron_interval = pybedtools.create_interval_from_list(
        ['chr1', '0', str(length), 'intron', '0', '+']
    )
    proxdist_dict = af.get_proxdist_from_intron(
        interval=intron_interval, distance=midpoint
    )
    assert 'prox' in proxdist_dict.keys()  # found a prox intron region.
    assert len(proxdist_dict['dist']) == 0  # found no dist intron
    assert len(proxdist_dict['prox']) == 1  # just found one prox intron
    assert proxdist_dict['prox'][0] == pybedtools.create_interval_from_list(
        ['chr1', '0', str(length), 'proxintron5', '0', '+']
    )

def test_split_prox_dist_2():
    print("Tests the core functionality of assigning proximal and distal "
          "intron spaces. Should return one intron based on the "
          "specified distance.")
    length = 11
    midpoint = 5
    intron_interval = pybedtools.create_interval_from_list(
        ['chr1', '0', str(length), 'intron', '0', '+']
    )
    proxdist_dict = af.get_proxdist_from_intron(
        interval=intron_interval, distance=midpoint
    )
    assert 'prox' in proxdist_dict.keys()  # found a prox intron region.
    assert 'dist' in proxdist_dict.keys()  # found a prox intron region.
    assert len(proxdist_dict['dist']) == 1  # found one dist intron
    assert len(proxdist_dict['prox']) == 2  # found two prox introns
    assert proxdist_dict['prox'][0] == pybedtools.create_interval_from_list(
        ['chr1', '0', '5', 'proxintron5', '0', '+']
    )
    assert proxdist_dict['prox'][1] == pybedtools.create_interval_from_list(
        ['chr1', '6', '11', 'proxintron5', '0', '+']
    )
    assert proxdist_dict['dist'][0] == pybedtools.create_interval_from_list(
        ['chr1', '5', '6', 'distintron5', '0', '+']
    )

def test_proxdist_regions_1():
    print("This tests that the a.create_definitions() function correctly "
          "identifies an intron (251-1250, 1-based inclusive; 250-1250, "
          "1-based inclusive) present in two genes, and that it "
          "correctly splits them into 1 proxintron and 2 proxintrons/1 dist"
          "intron respectively.")
    db_file = get_proxdist_db()
    chroms = ['chr1']
    species = 'hg19'
    append_chr = False
    fuzzy = False
    num_proxintrons = 0
    exons_dict, transcripts_dict, \
    cds_dict, features_dict, keys = a.create_definitions(
        db_file, chroms=chroms, species=species, append_chr=append_chr,
        fuzzy=fuzzy
    )
    for hashkey in features_dict.keys():
        for feature in features_dict[hashkey]:
            if feature.featuretype == 'distintron500':
                assert feature.attributes['gene_id'][0] == 'DISTINTRON_GENE'
            elif feature.featuretype == 'proxintron500':
                num_proxintrons += 1
    assert num_proxintrons == 3
