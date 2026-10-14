"""Exercise documented CLI entry points with small local input files."""

import gffutils
import pandas as pd
import pybedtools
import pytest

from annotator import annotator, build_gffutils_db, create_AS_STRUCTURE
from annotator import create_region_bedfiles, gene_name2id, miRNA_name2id


@pytest.fixture
def database(tmp_path):
    gtf = tmp_path / 'input.gtf'
    gtf.write_text(''.join('chr1\tx\t{}\t{}\t{}\t.\t+\t.\tgene_id "g1"; gene_name "Alpha"; transcript_id "tx1"; transcript_type "protein_coding"; gene_type "protein_coding";\n'.format(*r) for r in [
        ('gene',1,100), ('transcript',1,100), ('exon',1,100), ('CDS',20,80), ('UTR',1,19), ('UTR',81,100)]))
    path = tmp_path / 'input.db'
    build_gffutils_db.build_db(str(gtf), str(path), True, True, True)
    return gtf, path


def test_database_cli_builds_queryable_feature_database(database, tmp_path, monkeypatch):
    gtf, _ = database
    out = tmp_path / 'cli.db'
    monkeypatch.setattr('sys.argv', ['build_gffutils_db', '--annotation_file',str(gtf),'--db_file',str(out),
                                  '--disable_infer_genes','--disable_infer_transcripts'])
    build_gffutils_db.main()
    db = gffutils.FeatureDB(str(out))
    assert db.count_features_of_type('gene') == 1
    assert db.count_features_of_type('CDS') == 1
    assert [(f.start,f.end) for f in db.features_of_type('CDS')] == [(20,80)]


@pytest.mark.parametrize('by_transcript,name', [(False,'g1'),(True,'tx1')])
def test_region_cli_writes_only_requested_exports(database, tmp_path, monkeypatch, by_transcript, name):
    _, db = database
    cds = tmp_path / 'cds.bed'; exons = tmp_path / 'exons.bed'
    args = ['create_region_bedfiles','--db_file',str(db),'--species','gencode',
            '--cds_out',str(cds),'--exons_out',str(exons)]
    if by_transcript: args.append('--by-transcript')
    monkeypatch.setattr('sys.argv',args)
    create_region_bedfiles.main()
    assert [(x.start,x.end,x.name) for x in pybedtools.BedTool(str(cds))] == [(19,80,name)]
    assert [(x.start,x.end,x.name) for x in pybedtools.BedTool(str(exons))] == [(0,100,name)]


def test_annotation_cli_creates_coding_and_intergenic_rows(database, tmp_path, monkeypatch):
    _,db = database
    source = tmp_path / 'queries.bed'; out = tmp_path / 'annotated.bed'
    source.write_text('chr1\t30\t40\tcoding\t0\t+\nchr1\t200\t210\toutside\t0\t+\n')
    monkeypatch.setattr('sys.argv',['annotator','--input',str(source),'--output',str(out),'--gtfdb',str(db)])
    annotator.main()
    rows = [line.split('\t') for line in out.read_text().splitlines()]
    assert rows[0][:6] == ['chr1','30','40','coding','0','+']
    assert rows[0][6:9] == ['g1','Alpha','CDS']
    assert rows[1][6:9] == ['intergenic'] * 3


def test_gene_cli_converts_custom_mapping_and_preserves_measurement(tmp_path, monkeypatch):
    source = tmp_path / 'input.tsv'; source.write_text('symbol\tcount\nAlpha\t7\n')
    mapping = tmp_path / 'map.tsv'; mapping.write_text('g1\tAlpha\n')
    out = tmp_path / 'output.tsv'
    monkeypatch.setattr('sys.argv',['gene_name2id','--input',str(source),'--output',str(out),
                                  '--custom',str(mapping),'--name_col','symbol','--delete_original_column'])
    gene_name2id.main()
    assert pd.read_table(out).to_dict('list') == {'gene id':['g1'], 'count':[7]}


def test_mirna_cli_converts_custom_precursor_and_mature_ids(tmp_path, monkeypatch):
    source = tmp_path / 'input.tsv'; source.write_text('symbol\tcount\nmir-one\t7\n')
    mapping = tmp_path / 'map.tsv'; mapping.write_text('pre1\tmir-one\tmat1|mat2\n')
    out = tmp_path / 'output.tsv'
    monkeypatch.setattr('sys.argv',['miRNA_name2id','--input',str(source),'--output',str(out),
                                  '--custom',str(mapping),'--name_col','symbol','--add_mature'])
    miRNA_name2id.main()
    assert pd.read_table(out).to_dict('list') == {'miRNA id':['pre1'], 'count':[7], 'mature':['mat1|mat2']}


def test_as_structure_cli_writes_measured_gene_attributes(database, tmp_path, monkeypatch):
    _,db = database
    out = tmp_path / 'structure.gff'
    monkeypatch.setattr('sys.argv',['create_as_structure','--db_file',str(db),'--out_file',str(out),'--species','gencode'])
    create_AS_STRUCTURE.main()
    assert out.read_text().split('\t')[-1].strip() == 'gene_id=g1;mrna_length=100;premrna_length=100'


@pytest.mark.parametrize('module', [annotator, build_gffutils_db, create_region_bedfiles,
                                    create_AS_STRUCTURE, gene_name2id, miRNA_name2id])
def test_cli_help_lists_options_and_exits_successfully(module, monkeypatch, capsys):
    monkeypatch.setattr('sys.argv', ['tool', '--help'])
    with pytest.raises(SystemExit) as exc:
        module.main()
    assert exc.value.code == 0
    assert 'usage:' in capsys.readouterr().out


def test_annotation_cli_applies_prefix_strand_fallback_and_priority_files(tmp_path, monkeypatch):
    gtf=tmp_path/'input.gtf'
    gtf.write_text(''.join('1\tx\t{}\t{}\t{}\t.\t-\t.\tgene_id "g1"; gene_name "Alpha"; transcript_id "tx1"; transcript_type "protein_coding"; gene_type "protein_coding";\n'.format(*r) for r in [
        ('gene',1,100),('transcript',1,100),('exon',1,100),('CDS',20,80)]))
    db=tmp_path/'input.db';build_gffutils_db.build_db(str(gtf),str(db),True,True,True)
    source=tmp_path/'queries.bed';source.write_text('chr1\t30\t40\tquery\t0\t+\n')
    priority=tmp_path/'priority.txt';priority.write_text('protein_coding,CDS\nprotein_coding,exon\n')
    out=tmp_path/'output.bed'
    monkeypatch.setattr('sys.argv',['annotator','--input',str(source),'--output',str(out),'--gtfdb',str(db),
                                  '--append-chr','--unstranded','--limit-chroms-to','1',
                                  '--transcript-priority-file',str(priority),'--gene-priority-file',str(priority)])
    annotator.main()
    assert out.read_text().strip().split('\t')[6:9]==['g1','Alpha','CDS']


@pytest.mark.parametrize('strand,expected', [('+','miRNA'),('-','miRNA')])
def test_annotation_cli_identifies_mirna_exons(tmp_path,monkeypatch,strand,expected):
    gtf=tmp_path/'input.gtf'
    gtf.write_text(''.join('chr1\tx\t{}\t1\t100\t.\t{}\t.\tgene_id "g1"; gene_name "mir-one"; transcript_id "tx1"; transcript_type "miRNA"; gene_type "miRNA";\n'.format(kind,strand) for kind in ['gene','transcript','exon']))
    db=tmp_path/'input.db';build_gffutils_db.build_db(str(gtf),str(db),True,True,True)
    source=tmp_path/'queries.bed';source.write_text('chr1\t30\t40\tquery\t0\t'+strand+'\n')
    out=tmp_path/'output.bed'
    monkeypatch.setattr('sys.argv',['annotator','--input',str(source),'--output',str(out),'--gtfdb',str(db)])
    annotator.main()
    assert out.read_text().strip().split('\t')[6:9]==['g1','mir-one',expected]
