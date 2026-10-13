"""Offline GFF and tabular fixtures for gene and microRNA identifier conversion."""

import gffutils
import pandas as pd
import pytest

from annotator import gene_name2id as genes
from annotator import miRNA_name2id as mirnas


@pytest.fixture
def feature_db(tmp_path):
    path = tmp_path / 'features.db'
    gffutils.create_db('chr1\tx\tgene\t1\t50\t.\t+\t.\tID=g1;gene_id=g1;gene_name=Alpha\n'
                      'chr1\tx\tgene\t70\t90\t.\t+\t.\tID=g2;gene_id=g2;gene_name=Alpha\n'
                      'chr1\tx\tmiRNA_primary_transcript\t100\t160\t.\t+\t.\tID=pre1;Name=mir-one\n'
                      'chr1\tx\tmiRNA_primary_transcript\t170\t200\t.\t-\t.\tID=pre2;Name=mir-one\n'
                      'chr1\tx\tmiRNA\t110\t130\t.\t+\t.\tID=mat1;Derives_from=pre1\n'
                      'chr1\tx\tmiRNA\t140\t150\t.\t+\t.\tID=mat2;Derives_from=pre1\n',
                      dbfn=str(path), from_string=True)
    return str(path)


def test_gene_database_lookup_accumulates_ids_for_shared_name(feature_db):
    assert dict(genes.build_name_to_id_dict(feature_db, None)) == {'Alpha': ['g1', 'g2']}


def test_gene_custom_mapping_extends_database_and_ignores_comments(feature_db, tmp_path):
    custom = tmp_path / 'names.tsv'
    custom.write_text('# ID\tname\ng3\tAlpha\ng4\tBeta\n')
    assert dict(genes.build_name_to_id_dict(feature_db, str(custom))) == {
        'Alpha': ['g1', 'g2', 'g3'], 'Beta': ['g4']}


def test_precursor_database_lookup_accumulates_shared_names(feature_db):
    assert dict(mirnas.build_name_to_precursor_id_dict(feature_db, None)) == {'mir-one': ['pre1', 'pre2']}


def test_precursor_custom_mapping_uses_first_two_columns(tmp_path):
    custom = tmp_path / 'mirnas.tsv'
    custom.write_text('# comment\npre1\tmir-one\textra\npre2\tmir-one\n')
    assert dict(mirnas.build_name_to_precursor_id_dict(None, str(custom))) == {'mir-one': ['pre1', 'pre2']}


def test_mature_database_lookup_groups_ids_by_parent(feature_db):
    assert dict(mirnas.build_precursor_id_to_mature_id_dict(feature_db, None)) == {'pre1': ['mat1', 'mat2']}


def test_mature_custom_mapping_splits_pipe_delimited_ids(tmp_path):
    custom = tmp_path / 'mirnas.tsv'
    custom.write_text('# comment\npre1\tmir-one\tmat1|mat2\n')
    assert dict(mirnas.build_precursor_id_to_mature_id_dict(None, str(custom))) == {'pre1': ['mat1', 'mat2']}


@pytest.mark.parametrize('module,function', [(genes, 'build_name_to_id_dict'),
                                            (mirnas, 'build_name_to_precursor_id_dict'),
                                            (mirnas, 'build_precursor_id_to_mature_id_dict')])
def test_identifier_mapping_without_sources_is_empty(module, function):
    assert getattr(module, function)(None, None) == {}


@pytest.mark.parametrize('module,function', [(genes, 'get_id_from_namerow'),
                                            (mirnas, 'get_precursor_id_from_namerow'),
                                            (mirnas, 'get_mature_id_from_idrow')])
def test_row_lookup_joins_ids_and_handles_unknown_name(module, function):
    lookup = getattr(module, function)
    assert lookup({'symbol': 'known'}, {'known': ['a', 'b']}, 'symbol') == 'a|b'
    assert lookup({'symbol': 'missing'}, {'known': ['a']}, 'symbol') == ''


def test_gene_row_lookup_normalizes_numeric_names_to_strings():
    assert genes.get_id_from_namerow({'symbol': 123}, {'123': ['g1']}, 'symbol') == 'g1'


@pytest.mark.parametrize('module,function,column', [(genes, 'add_id_to_dataframe', 'gene id'),
                                                   (mirnas, 'add_precursor_id_to_dataframe', 'miRNA id'),
                                                   (mirnas, 'add_mature_id_to_dataframe', 'mature')])
def test_dataframe_lookup_preserves_rows_and_unrelated_columns(module, function, column):
    df = pd.DataFrame({'symbol': ['known', 'missing'], 'count': [7, 2]})
    assert getattr(module, function)(df, {'known': ['a', 'b']}, 'symbol') is df
    assert df.to_dict('list') == {'symbol': ['known', 'missing'], 'count': [7, 2], column: ['a|b', '']}


@pytest.mark.parametrize('module,function,column', [(genes, 'add_id_to_dataframe', 'gene id'),
                                                   (mirnas, 'add_precursor_id_to_dataframe', 'miRNA id')])
def test_dataframe_lookup_refuses_to_replace_existing_id_column(module, function, column):
    df = pd.DataFrame({'symbol': ['known'], column: ['original']})
    with pytest.raises(AssertionError):
        getattr(module, function)(df, {'known': ['new']}, 'symbol')
    assert df[column].tolist() == ['original']


@pytest.mark.parametrize('delete_original', [False, True])
def test_gene_table_conversion_preserves_measurements_and_selects_index(tmp_path, delete_original):
    source = tmp_path / 'counts.tsv'; source.write_text('symbol\tcount\nAlpha\t7\nBeta\t2\n')
    mapping = tmp_path / 'names.tsv'; mapping.write_text('g1\tAlpha\ng2\tBeta\n')
    df = genes.read_and_append_id(str(source), '\t', 'symbol', None, str(mapping), delete_original)
    assert df.index.name == ('gene id' if delete_original else 'symbol')
    assert df.index.tolist() == (['g1', 'g2'] if delete_original else ['Alpha', 'Beta'])
    assert df['count'].tolist() == [7, 2]
    out = tmp_path / 'out.tsv'
    genes.convert(str(source), '\t', 'symbol', str(out), None, str(mapping), delete_original)
    pd.testing.assert_frame_equal(pd.read_csv(out, sep='\t', index_col=0), df)


@pytest.mark.parametrize('delete_original', [False, True])
def test_mirna_table_conversion_can_include_mature_ids(tmp_path, delete_original):
    source = tmp_path / 'counts.tsv'; source.write_text('symbol\tcount\nmir-one\t7\nmir-two\t2\n')
    mapping = tmp_path / 'names.tsv'; mapping.write_text('pre1\tmir-one\tmat1|mat2\npre2\tmir-two\tmat3\n')
    df = mirnas.read_and_append_id(str(source), '\t', 'symbol', None, str(mapping), True, delete_original)
    assert df.index.name == ('miRNA id' if delete_original else 'symbol')
    assert df.index.tolist() == (['pre1', 'pre2'] if delete_original else ['mir-one', 'mir-two'])
    assert df['mature'].tolist() == ['mat1|mat2', 'mat3']
    assert df['count'].tolist() == [7, 2]
    out = tmp_path / 'out.tsv'
    mirnas.convert(str(source), '\t', 'symbol', str(out), None, str(mapping), True, delete_original)
    pd.testing.assert_frame_equal(pd.read_csv(out, sep='\t', index_col=0), df)


@pytest.mark.parametrize('module,function', [(genes, 'build_name_to_id_dict'),
                                            (mirnas, 'build_name_to_precursor_id_dict'),
                                            (mirnas, 'build_precursor_id_to_mature_id_dict')])
def test_mapping_sources_report_missing_paths_without_losing_available_custom_data(module, function, tmp_path, capsys):
    missing = tmp_path / 'missing.db'
    assert getattr(module, function)(str(missing), str(tmp_path / 'missing.tsv')) == {}
    out = capsys.readouterr().out
    assert 'Could not use the gffdb file' in out
    assert 'Could not use the custom file' in out
