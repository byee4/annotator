"""Small synthetic examples for annotation coordinates and format contracts."""

from collections import defaultdict
from types import SimpleNamespace

import gffutils
import pybedtools
import pytest

from annotator import annotation_functions as af
from annotator import annotate_bed as ab


def feature(start, end, strand='+', transcripts=('tx',), featuretype='CDS'):
    return gffutils.Feature(seqid='chr1', start=start, end=end, strand=strand,
                            featuretype=featuretype,
                            attributes={'transcript_id': list(transcripts)})


@pytest.mark.parametrize('source,expected', [
    ('ensembl', ('CDS', None, 'three_prime_utr', 'five_prime_utr', 'transcript_biotype', 'gene_biotype')),
    ('ce10', ('CDS', 'UTR', None, None, 'transcript_biotype', 'gene_biotype')),
    ('wormbase', ('CDS', None, 'three_prime_utr', 'five_prime_utr', 'transcript_biotype', 'gene_biotype')),
    ('ce11', ('CDS', None, 'three_prime_utr', 'five_prime_utr', 'transcript_biotype', 'gene_biotype')),
    ('gencode', ('CDS', 'UTR', None, None, 'transcript_type', 'gene_type')),
    ('unknown', ('CDS', 'UTR', None, None, 'transcript_type', 'gene_type')),
    ('genbank', ('Genbank CDS', None, None, None, 'gene_biotype', 'gene_biotype')),
    ('dm6_genbank', ('Genbank CDS', None, None, None, 'gene_biotype', 'gene_biotype')),
    ('refseq', ('CDS', None, None, None, 'gene_biotype', 'gene_biotype')),
    ('ucsc', ('CDS', None, '3UTR', '5UTR', 'gene_biotype', 'gene_type')),
])
def test_get_keys_selects_source_vocabulary(source, expected):
    keys = af.get_keys(source)
    assert tuple(keys[k] for k in ('cds', 'utr', 'utr3', 'utr5', 'transcript_type', 'gene_type')) == expected
    assert set(keys) == {'cds', 'utr', 'utr3', 'utr5', 'transcript_type', 'gene_type',
                         'gene', 'gene_name', 'transcript', 'transcript_id', 'exon', 'gene_id'}


def test_get_keys_returns_independent_mappings():
    first = af.get_keys('gencode')
    first['cds'] = 'changed'
    assert af.get_keys('gencode')['cds'] == 'CDS'


def test_get_all_cds_dict_unions_extrema_for_each_transcript():
    db = SimpleNamespace(features_of_type=lambda kind: iter([
        feature(100, 110, transcripts=('tx', 'shared')),
        feature(50, 75), feature(105, 150), feature(80, 90),
    ]) if kind == 'coding' else iter([]))
    assert af.get_all_cds_dict(db, 'coding') == {
        'tx': {'low': 50, 'hi': 150}, 'shared': {'low': 100, 'hi': 110}}


def test_get_all_cds_dict_empty_database():
    db = SimpleNamespace(features_of_type=lambda kind: iter([]))
    assert af.get_all_cds_dict(db, 'CDS') == {}


def test_chromosome_set_deduplicates_database_rows():
    db = gffutils.create_db('chr1\tx\tgene\t1\t3\t.\t+\t.\tgene_id "g1";\n'
                           'chr2\tx\tgene\t1\t3\t.\t-\t.\tgene_id "g2";\n'
                           'chr1\tx\tgene\t8\t9\t.\t+\t.\tgene_id "g3";\n',
                           dbfn=':memory:', from_string=True)
    assert af.chromosome_set(db) == {'chr1', 'chr2'}


@pytest.mark.parametrize('strand,start,end,expected', [
    ('+', 10, 49, '5utr'), ('+', 101, 120, '3utr'),
    ('-', 10, 49, '3utr'), ('-', 101, 120, '5utr'),
    ('+', 60, 80, 'unclassified_utr'), ('.', 10, 49, 'unclassified_utr'),
])
def test_classify_utr_respects_coding_span_and_strand(strand, start, end, expected):
    assert af.classify_utr(feature(start, end, strand), {'tx': {'low': 50, 'hi': 100}}) == expected


def test_classify_utr_without_coding_transcript_is_unclassified():
    assert af.classify_utr(feature(1, 10), {}) == 'unclassified_utr'


@pytest.mark.parametrize('length,distance,prox,dist', [
    (100, 50, [(10, 110)], []),
    (101, 50, [(10, 60), (61, 111)], [(60, 61)]),
    (200, 50, [(10, 60), (160, 210)], [(60, 160)]),
])
def test_proximal_distal_split_preserves_half_open_coordinates(length, distance, prox, dist):
    intron = pybedtools.create_interval_from_list(['chr2', '10', str(10 + length), 'tx', '7', '-'])
    result = af.get_proxdist_from_intron(intron, distance)
    assert [(x.start, x.end) for x in result['prox']] == prox
    assert [(x.start, x.end) for x in result['dist']] == dist
    assert sum(len(x) for values in result.values() for x in values) == length
    assert all((x.chrom, x.strand, x.score) == ('chr2', '-', '7') for values in result.values() for x in values)


@pytest.mark.parametrize('exons,expected', [
    ([{'start': 1, 'end': 20}, {'start': 31, 'end': 50}], [{'start': 21, 'end': 30}]),
    ([{'start': 31, 'end': 50}, {'start': 1, 'end': 20}], [{'start': 21, 'end': 30}]),
    ([{'start': 1, 'end': 50}], []),
    ([], []),
])
def test_find_introns_uses_inclusive_exon_coordinates(exons, expected):
    assert af.find_introns({'start': 1, 'end': 50}, exons) == expected


@pytest.mark.parametrize('tx_type,gene_type,distinct,expected', [
    ('protein_coding', 'lncRNA', None, 'protein_coding'),
    ('miRNA', 'protein_coding', None, 'miRNA'),
    ('retained_intron', 'protein_coding', None, 'non_coding'),
    ('lncRNA', 'lncRNA', ['lncRNA'], 'lncRNA'),
    ('lncRNA', 'lncRNA', ['snoRNA'], 'non_coding'),
])
def test_classify_transcript_type_preserves_selected_noncoding_classes(tx_type, gene_type, distinct, expected):
    assert ab.classify_transcript_type(tx_type, gene_type, distinct) == expected


@pytest.mark.parametrize('value,expected', [('protein_coding', True), ('miRNA', False), ('', False)])
def test_is_protein_coding_requires_exact_type(value, expected):
    assert ab.is_protein_coding(value) is expected


def test_split_bed_interval_converts_only_coordinates():
    assert ab.split_bed_interval('chr2\t0\t10\tpeak\t3.5\t-\n') == ('chr2', 0, 10, 'peak', '3.5', '-')


@pytest.mark.parametrize('line,delimiter,expected', [
    ('intergenic', ',', ('intergenic',) * 3),
    ('tx,10,20,+,CDS,g1,Alpha,protein_coding,protein_coding,overlap', ',', ('CDS', 'g1', 'Alpha')),
    ('tx:10:20:+:exon:g2:Beta:type:type:overlap', ':', ('exon', 'g2', 'Beta')),
])
def test_split_string_extracts_public_annotation_columns(line, delimiter, expected):
    assert ab.split_string(line, delimiter) == expected


@pytest.mark.parametrize('text,expected', [('a|b', ['a', 'b']), ('a|', ['a', '']), ('', [''])])
def test_parse_annotation_string_preserves_feature_order(text, expected):
    assert ab.parse_annotation_string(text) == expected


@pytest.mark.parametrize('strand,query,expected', [
    ('+', (9, 20), ['same']), ('-', (9, 20), ['other']),
    ('+', (0, 10), ['same']), ('+', (20, 30), []),
])
def test_stranded_overlap_observes_bed_to_gtf_boundary(strand, query, expected):
    features = defaultdict(list)
    a = feature(10, 20, '+'); a.id = 'same'
    b = feature(10, 20, '-'); b.id = 'other'
    features['chr1', 0, '+'] = [a]; features['chr1', 0, '-'] = [b]
    result = ab.get_all_overlapping_features_from_query_stranded('chr1', *query, strand, features)
    assert [x.id for x in result] == expected


@pytest.mark.parametrize('strand,features_on_plus,expected', [
    (None, True, ['plus', 'minus']), ('+', True, ['plus']), ('+', False, ['minus']),
])
def test_unstranded_queries_prefer_requested_strand_then_fall_back(strand, features_on_plus, expected):
    features = defaultdict(list)
    a = feature(10, 20, '+'); a.id = 'plus'
    b = feature(10, 20, '-'); b.id = 'minus'
    if features_on_plus: features['chr1', 0, '+'] = [a]
    features['chr1', 0, '-'] = [b]
    result = ab.get_all_overlapping_features_from_query('chr1', 9, 20, strand, features, stranded=False)
    assert [x.id for x in result] == expected


@pytest.mark.parametrize('delimiter,text', [(',', 'protein_coding,CDS\nmiRNA,exon\n'),
                                            ('\t', 'protein_coding\tCDS\nmiRNA\texon\n')])
def test_priority_file_preserves_class_and_region_order(tmp_path, delimiter, text):
    from annotator.annotator import parse_annotation_priority
    path = tmp_path / 'priority.txt'
    path.write_text(text)
    assert parse_annotation_priority(str(path), delimiter) == [['protein_coding','CDS'], ['miRNA','exon']]


def test_profiling_decorator_preserves_arguments_result_and_prints_stats(capsys):
    @ab.do_cprofile
    def add(left, *, right):
        return left + right
    assert add(2, right=3) == 5
    assert 'function calls' in capsys.readouterr().out


def test_profiling_decorator_propagates_failure_and_prints_stats(capsys):
    @ab.do_cprofile
    def failing():
        raise ValueError('sentinel')
    with pytest.raises(ValueError, match='sentinel'):
        failing()
    assert 'function calls' in capsys.readouterr().out
