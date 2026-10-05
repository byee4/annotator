"""Keep the committed Ensembl example database in sync with its source GTF."""

from pathlib import Path

import gffutils

from annotator import annotate_bed
from annotator.annotator import GENE_PRIORITY, TRANSCRIPT_PRIORITY
from annotator.build_gffutils_db import build_db


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "ensembl_gtf"
GTF = EXAMPLE / "mini_ensembl.gtf"
DB = EXAMPLE / "mini_ensembl.gtf.db"
BED = EXAMPLE / "queries.bed"


def feature_records(db):
    return sorted(
        (
            feature.featuretype,
            feature.seqid,
            feature.start,
            feature.end,
            feature.strand,
            tuple(sorted((key, tuple(value)) for key, value in feature.attributes.items())),
        )
        for feature in db.all_features()
    )


def test_committed_ensembl_database_matches_gtf(tmp_path):
    rebuilt_file = tmp_path / "rebuilt.gtf.db"
    build_db(str(GTF), str(rebuilt_file))

    committed = gffutils.FeatureDB(str(DB))
    rebuilt = gffutils.FeatureDB(str(rebuilt_file))
    assert feature_records(committed) == feature_records(rebuilt)
    assert len(feature_records(committed)) == 7
    assert list(committed.features_of_type("gene")) == []


def test_ensembl_example_annotates_coding_utr_and_intergenic_queries(tmp_path):
    output = tmp_path / "annotated.bed"
    annotate_bed.annotate_bed(
        str(DB), [str(BED)], [str(output)], True, [],
        [row[:] for row in TRANSCRIPT_PRIORITY],
        [row[:] for row in GENE_PRIORITY],
        "ensembl", False, 0, 1,
    )

    rows = [line.rstrip("\n").split("\t") for line in output.read_text().splitlines()]
    assert [(row[3], row[6], row[8]) for row in rows] == [
        ("example_cds", "ENSG_EXAMPLE_1", "CDS"),
        ("example_5utr", "ENSG_EXAMPLE_1", "five_prime_utr"),
        ("example_intergenic", "intergenic", "intergenic"),
        ("issue_3_transcript", "ENSG00000272512", "noncoding_exon"),
    ]
