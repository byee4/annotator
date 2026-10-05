"""Regression tests for coding and noncoding GTF biotype classification."""

import gffutils
import pytest

from annotator import annotate_bed
from annotator.annotator import GENE_PRIORITY, TRANSCRIPT_PRIORITY


def annotate_features(tmp_path, gtf_format, gene_key, transcript_key, transcript_type,
                      include_cds=True):
    gene_attrs = (
        f'gene_id "ENSG1"; gene_name "Example"; '
        f'{gene_key} "protein_coding";'
    )
    transcript_attrs = (
        f'gene_id "ENSG1"; transcript_id "ENST1"; gene_name "Example"; '
        f'{gene_key} "protein_coding"; {transcript_key} "{transcript_type}";'
    )
    features = [("gene", gene_attrs), ("transcript", transcript_attrs),
                ("exon", transcript_attrs)]
    if include_cds:
        features.append(("CDS", transcript_attrs))
    gtf = tmp_path / "reference.gtf"
    gtf.write_text("\n".join(
        f"chr1\ttest\t{feature}\t100\t200\t.\t+\t.\t{attributes}"
        for feature, attributes in features
    ) + "\n")
    db_file = tmp_path / "reference.gtf.db"
    gffutils.create_db(
        str(gtf), dbfn=str(db_file), force=True,
        keep_order=True, merge_strategy="merge", sort_attribute_values=True,
        disable_infer_genes=False, disable_infer_transcripts=False,
    )
    _, _, cds, indexed_features, keys = annotate_bed.create_definitions(
        str(db_file), gtf_format=gtf_format,
    )
    return annotate_bed.annotate(
        "chr1", 119, 120, "query", "0", "+", True,
        [priority[:] for priority in TRANSCRIPT_PRIORITY],
        [priority[:] for priority in GENE_PRIORITY],
        indexed_features, cds, keys,
    )


@pytest.mark.parametrize(
    "gtf_format,gene_key,transcript_key",
    [
        pytest.param("gencode", "gene_type", "transcript_type",
                     id="gencode-standard-keys"),
        pytest.param("gencode", "gene_biotype", "transcript_biotype",
                     id="gencode-setting-with-biotype-keys"),
        pytest.param("ensembl", "gene_type", "transcript_type",
                     id="ensembl-setting-with-type-keys"),
    ],
)
def test_protein_coding_cds_is_selected(
    tmp_path, gtf_format, gene_key, transcript_key,
):
    result = annotate_features(
        tmp_path, gtf_format, gene_key, transcript_key, "protein_coding",
    )
    assert result[6:9] == ("ENSG1", "Example", "CDS")


def test_noncoding_transcript_in_coding_gene_stays_noncoding(tmp_path):
    result = annotate_features(
        tmp_path, "gencode", "gene_type", "transcript_type",
        "retained_intron", include_cds=False,
    )
    assert result[8] == "noncoding_exon"
