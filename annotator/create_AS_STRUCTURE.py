#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse
import gffutils
from tqdm import trange
from . import annotation_functions as af

def create_as_structure(db, species, out_file):
    """
    
    :param db: gffutils.FeatureDB
    :param species: string
        sets the gtf nomenclature
    :param out_file: string
    :return: 
    """
    keys = af.get_keys(species)
    genes = af.get_gene_to_transcript_dict(
        db=db,
        gene_id_key=keys['gene_id'],
        transcript_id_key=keys['transcript_id']
    )
    exons = af.get_all_exons_dict(
        db=db,
        exon_key=keys['exon'],
        transcript_id_key=keys['transcript_id']
    )
    transcripts = af.get_all_transcripts_dict(
        db=db,
        transcript_key=keys['transcript'],
        transcript_id_key=keys['transcript_id']
    )
    upstream_downstream = af.most_upstream_downstream_positions(
        genes_dict=genes,
        transcripts_dict=transcripts
    )
    progress = trange(db.count_features_of_type(keys['gene']))
    with open(out_file, 'w') as f:
        for gene in db.features_of_type(keys['gene']):
            gene_id = gene.attributes[keys['gene_id']][0]  # just take first gene_id
            f.write(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gene.seqid,
                    "AS_STRUCTURE",
                    "gene",
                    upstream_downstream[gene_id]['start'],
                    upstream_downstream[gene_id]['end'],
                    '.',
                    gene.strand,
                    '.',
                    "gene_id={};mrna_length={};premrna_length={}".format(
                        gene_id,
                        af.get_mrna_lengths(
                            gene_id=gene_id,
                            exons_dict=exons,
                            genes_dict=genes
                        ),
                        af.get_premrna_lengths(
                            gene_id=gene_id,
                            upstream_downstream=upstream_downstream,
                        )
                    )
                )
            )
            progress.update(1)


def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db_file",
        required=True,
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )
    parser.add_argument(
        "--species",
        required=False,
        default='hg19'
    )
    args = parser.parse_args()

    out_file = args.out_file
    db_file = args.db_file
    species = args.species

    db = gffutils.FeatureDB(db_file)
    create_as_structure(db=db, species=species, out_file=out_file)

if __name__ == "__main__":
    main()
