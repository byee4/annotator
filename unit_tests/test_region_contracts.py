"""Small real BED/GTF fixtures for region exports and transcript geometry."""

import gffutils
import pybedtools
import pytest

from annotator import annotation_functions as af
from annotator import create_region_bedfiles as regions
from annotator import create_AS_STRUCTURE as structure


@pytest.fixture
def db():
    rows = [('gene', 1, 1600), ('transcript', 1, 1600),
            ('exon', 1, 100), ('exon', 1501, 1600),
            ('CDS', 20, 100), ('CDS', 1501, 1550),
            ('UTR', 1, 19), ('UTR', 1551, 1600)]
    text = ''.join('chr1\tx\t{}\t{}\t{}\t.\t+\t.\tgene_id "g1"; gene_name "Alpha"; transcript_id "tx1";\n'.format(*r) for r in rows)
    return gffutils.create_db(text, dbfn=':memory:', from_string=True,
                             disable_infer_genes=True, disable_infer_transcripts=True)


@pytest.fixture
def keys():
    return af.get_keys('gencode')


def coords(tool):
    return [(x.chrom, x.start, x.end, x.name, x.strand) for x in tool]


@pytest.mark.parametrize('by_transcript,name', [(False, 'g1'), (True, 'tx1')])
def test_create_bedtools_converts_inclusive_gtf_to_half_open_bed(db, keys, by_transcript, name):
    assert coords(regions.create_bedtools(list(db.features_of_type('exon')), keys, by_transcript)) == [
        ('chr1', 0, 100, name, '+'), ('chr1', 1500, 1600, name, '+')]


def test_create_bedtools_retains_each_associated_id(keys):
    f = gffutils.Feature(seqid='chr2', start=5, end=5, strand='-',
                         attributes={'gene_id': ['a', 'b']})
    assert coords(regions.create_bedtools([f], keys)) == [
        ('chr2', 4, 5, 'a', '-'), ('chr2', 4, 5, 'b', '-')]


def test_merge_bedtool_by_gene_preserves_gene_and_strand_boundaries():
    bed = pybedtools.BedTool('chr1\t10\t20\ta\t0\t+\nchr1\t15\t30\ta\t0\t+\n'
                            'chr1\t15\t40\tb\t0\t+\nchr1\t18\t25\ta\t0\t-\n'
                            'chr2\t10\t20\ta\t0\t+\n', from_string=True)
    df = regions.merge_bedtool_by_gene(bed)
    actual = set(df[['chrom', 'start', 'end', 'name', 'strand']].itertuples(index=False, name=None))
    assert actual == {('chr1', 10, 30, 'a', '+'), ('chr1', 15, 40, 'b', '+'),
                      ('chr1', 18, 25, 'a', '-'), ('chr2', 10, 20, 'a', '+')}


@pytest.mark.parametrize('by_transcript,name', [(False, 'g1'), (True, 'tx1')])
def test_cds_region_export_and_returned_bed_have_same_coordinates(db, keys, tmp_path, by_transcript, name):
    out = tmp_path / 'cds.bed'
    result = regions.create_cds_region_bedfile(db, keys, str(out), by_transcript)
    expected = [('chr1', 19, 100, name, '+'), ('chr1', 1500, 1550, name, '+')]
    assert coords(result) == expected
    assert coords(pybedtools.BedTool(str(out))) == expected


@pytest.mark.parametrize('by_transcript,name', [(False, 'g1'), (True, 'tx1')])
def test_exon_region_export_preserves_span_and_selected_name(db, keys, tmp_path, by_transcript, name):
    out = tmp_path / 'exon.bed'
    assert regions.create_exon_region_bedfile(db, keys, str(out), by_transcript) is None
    assert coords(pybedtools.BedTool(str(out))) == [('chr1', 0, 100, name, '+'), ('chr1', 1500, 1600, name, '+')]


def test_gene_region_export_uses_gene_identifier(db, keys, tmp_path):
    out = tmp_path / 'gene.bed'
    regions.create_gene_region_bedfile(db, keys, str(out))
    assert coords(pybedtools.BedTool(str(out))) == [('chr1', 0, 1600, 'g1', '+')]


@pytest.mark.parametrize('by_transcript,name', [(False, 'g1'), (True, 'tx1')])
def test_utr_exports_classify_generic_utr_features(db, keys, tmp_path, by_transcript, name):
    five = tmp_path / '5utr.bed'; three = tmp_path / '3utr.bed'
    result5, result3 = regions.create_utr_region_bedfiles(
        db, keys, af.get_all_cds_dict(db, 'CDS'), str(three), str(five), by_transcript)
    assert coords(result5) == [('chr1', 0, 19, name, '+')]
    assert coords(result3) == [('chr1', 1550, 1600, name, '+')]
    assert coords(pybedtools.BedTool(str(five))) == coords(result5)
    assert coords(pybedtools.BedTool(str(three))) == coords(result3)


def test_utr_exports_accept_explicit_ensembl_utr_types(keys):
    f5 = gffutils.Feature(seqid='chr2', start=1, end=10, strand='-', attributes={'transcript_id':['tx'], 'gene_id':['g']})
    f3 = gffutils.Feature(seqid='chr2', start=80, end=100, strand='-', attributes={'transcript_id':['tx'], 'gene_id':['g']})
    class DB:
        def features_of_type(self, kind):
            return iter({'five_prime_utr': [f5], 'three_prime_utr': [f3]}.get(kind, []))
    five, three = regions.create_utr_region_bedfiles(DB(), af.get_keys('ensembl'), {}, None, None, True)
    assert coords(five) == [('chr2', 0, 10, 'tx', '-')]
    assert coords(three) == [('chr2', 79, 100, 'tx', '-')]


def test_utr_exports_report_unsupported_vocabulary(capsys):
    assert regions.create_utr_region_bedfiles(None, af.get_keys('refseq'), {}, None, None, True) == 1
    assert 'Having trouble guessing' in capsys.readouterr().out


@pytest.mark.parametrize('by_transcript,name', [(False, 'g1'), (True, 'tx1')])
def test_intron_return_values_partition_inclusive_exon_gap(db, keys, by_transcript, name):
    exons = af.get_all_exons_dict(db, 'exon', 'transcript_id')
    transcripts = af.get_all_transcripts_dict(db, 'transcript', 'transcript_id')
    prox, dist, all_introns = regions.create_intron_region_bedfiles(
        db, exons, transcripts, keys, None, None, None, by_transcript)
    assert coords(prox) == [('chr1', 100, 600, name, '+'), ('chr1', 1000, 1500, name, '+')]
    assert coords(dist) == [('chr1', 600, 1000, name, '+')]
    assert coords(all_introns) == [('chr1', 100, 1500, name, '+')]


def test_gene_intron_files_match_returned_regions(db, keys, tmp_path):
    paths = [tmp_path / n for n in ['prox.bed', 'dist.bed', 'all.bed']]
    tools = regions.create_intron_region_bedfiles(
        db, af.get_all_exons_dict(db, 'exon', 'transcript_id'),
        af.get_all_transcripts_dict(db, 'transcript', 'transcript_id'), keys,
        *(str(p) for p in paths), False)
    assert [coords(pybedtools.BedTool(str(p))) for p in paths] == [coords(t) for t in tools]


def test_transcript_geometry_chooses_first_longest_and_outer_bounds():
    genes = {'g': ['short', 'long', 'tied']}
    transcripts = {'short': {'start': 10, 'end': 20}, 'long': {'start': 5, 'end': 50},
                   'tied': {'start': 7, 'end': 52}}
    assert af.get_longest_transcripts(genes, transcripts) == {'g': 'long'}
    outer = af.most_upstream_downstream_positions(genes, transcripts)
    assert outer == {'g': {'start': 5, 'end': 52}}
    assert af.get_premrna_lengths('g', outer) == 48


def test_mrna_length_counts_exon_union_across_isoforms():
    exon = lambda start,end: {'chrom':'chr1', 'start': start, 'end':end, 'strand':'+'}
    assert af.get_mrna_lengths('g', {'tx1':[exon(1,10),exon(31,40)],
                                    'tx2':[exon(5,15),exon(31,40)]}, {'g':['tx1','tx2']}) == 25


def test_as_structure_export_includes_mrna_and_premrna_lengths(db, tmp_path):
    out = tmp_path / 'structure.gff'
    structure.create_as_structure(db, 'gencode', str(out))
    fields = out.read_text().strip().split('\t')
    assert fields[:8] == ['chr1','AS_STRUCTURE','gene','1','1600','.','+','.']
    assert fields[8] == 'gene_id=g1;mrna_length=200;premrna_length=1600'


def test_total_cds_length_merges_overlapping_coding_bases(db, keys, capsys):
    from annotator import get_region_lengths as lengths
    assert lengths.calculate_total_cds_length(db, keys) == 131
    assert 'cds lengths (total): 131' in capsys.readouterr().out


def test_total_intron_length_reports_each_partition(db, keys, capsys):
    from annotator import get_region_lengths as lengths
    lengths.calculate_total_intron_lengths(
        db, af.get_all_exons_dict(db, 'exon', 'transcript_id'),
        af.get_all_transcripts_dict(db, 'transcript', 'transcript_id'), keys)
    out = capsys.readouterr().out
    assert 'proximal intron lengths (total): 1000' in out
    assert 'distal intron lengths (total): 400' in out
    assert 'all intron lengths (total): 1400' in out


def test_total_lengths_wrapper_reports_actual_coding_utr_and_intron_bases(db, tmp_path, capsys):
    from annotator import get_region_lengths as lengths
    path = tmp_path / 'database.db'
    db.conn.backup(__import__('sqlite3').connect(str(path)))
    assert lengths.calculate_total_lengths(str(path), 'gencode') is None
    out = capsys.readouterr().out
    for text in ['cds lengths (total): 131', 'five_prime_utr_lengths (total): 19',
                 'three_prime_utr_lengths (total): 50', 'all intron lengths (total): 1400']:
        assert text in out


def test_gene_export_preserves_distinct_features_with_repeated_gene_id(keys,tmp_path,capsys):
    genes=[gffutils.Feature(seqid='chr1',start=1,end=10,strand='+',attributes={'gene_id':['same']}),
           gffutils.Feature(seqid='chr2',start=20,end=30,strand='-',attributes={'gene_id':['same']})]
    class DB:
        def features_of_type(self,kind):
            assert kind=='gene'
            return iter(genes)
    out=tmp_path/'genes.bed'
    regions.create_gene_region_bedfile(DB(),keys,str(out))
    assert coords(pybedtools.BedTool(str(out)))==[('chr1',0,10,'same','+'),('chr2',19,30,'same','-')]
    assert 'hashkey exists for gene same' in capsys.readouterr().out
