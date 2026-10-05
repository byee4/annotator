"""Unit tests for the public annotation lookup helpers."""

from types import SimpleNamespace

from annotator.annotation_functions import (
    gene_id_to_name,
    get_all_exons_dict,
    get_all_transcripts_dict,
    get_gene_to_transcript_dict,
)


class FeatureDBStub:
    def __init__(self, features):
        self.features = features

    def features_of_type(self, feature_type):
        return iter(self.features.get(feature_type, ()))


def feature(attributes, start=1, end=10, seqid="chr1", strand="+"):
    return SimpleNamespace(
        attributes=attributes,
        start=start,
        end=end,
        seqid=seqid,
        strand=strand,
    )


def test_gene_to_transcript_dict_groups_transcripts_by_gene():
    db = FeatureDBStub(
        {
            "transcript": [
                feature({"gene_id": ["gene-a"], "transcript_id": ["tx-1"]}),
                feature({"gene_id": ["gene-a"], "transcript_id": ["tx-2"]}),
                feature({"gene_id": ["gene-b"], "transcript_id": ["tx-3"]}),
            ],
            "gene": [feature({"gene_id": ["unrelated"]})],
        }
    )

    assert get_gene_to_transcript_dict(db, "gene_id", "transcript_id") == {
        "gene-a": ["tx-1", "tx-2"],
        "gene-b": ["tx-3"],
    }


def test_gene_to_transcript_dict_uses_supplied_attribute_keys():
    db = FeatureDBStub(
        {
            "transcript": [
                feature({"gene": ["gene-a", "gene-b"], "isoform": ["tx-1", "tx-2"]})
            ]
        }
    )

    assert get_gene_to_transcript_dict(db, "gene", "isoform") == {
        "gene-a": ["tx-1", "tx-2"],
        "gene-b": ["tx-1", "tx-2"],
    }


def test_gene_to_transcript_dict_returns_empty_mapping_without_transcripts():
    db = FeatureDBStub({"gene": [feature({"gene_id": ["gene-a"]})]})

    assert get_gene_to_transcript_dict(db, "gene_id", "transcript_id") == {}


def test_all_transcripts_dict_preserves_coordinates_for_each_transcript_id():
    db = FeatureDBStub(
        {
            "mRNA": [
                feature({"ID": ["tx-1", "tx-2"]}, start=10, end=30),
                feature({"ID": ["tx-3"]}, start=40, end=70),
            ],
            "transcript": [feature({"ID": ["unrelated"]})],
        }
    )

    assert get_all_transcripts_dict(db, "mRNA", "ID") == {
        "tx-1": {"start": 10, "end": 30},
        "tx-2": {"start": 10, "end": 30},
        "tx-3": {"start": 40, "end": 70},
    }


def test_all_exons_dict_keeps_each_exon_and_its_strand():
    db = FeatureDBStub(
        {
            "exon": [
                feature({"transcript_id": ["tx-1", "tx-2"]}, 10, 20, "chr2", "-"),
                feature({"transcript_id": ["tx-1"]}, 30, 40, "chr2", "-"),
            ]
        }
    )

    assert get_all_exons_dict(db, "exon", "transcript_id") == {
        "tx-1": [
            {"chrom": "chr2", "start": 10, "end": 20, "strand": "-"},
            {"chrom": "chr2", "start": 30, "end": 40, "strand": "-"},
        ],
        "tx-2": [{"chrom": "chr2", "start": 10, "end": 20, "strand": "-"}],
    }


def test_gene_id_to_name_accepts_list_and_scalar_gene_ids():
    db = FeatureDBStub(
        {
            "gene": [
                feature({"gene_id": ["gene-a"], "symbol": ["Alpha"]}),
                feature({"gene_id": "gene-b", "symbol": ["Beta"]}),
            ],
            "transcript": [feature({"gene_id": ["unrelated"], "symbol": ["Other"]})],
        }
    )

    assert gene_id_to_name(db, "symbol") == {
        "gene-a": "Alpha",
        "gene-b": "Beta",
    }
