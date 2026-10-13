#!/usr/bin/env python

# Parses a GTF file and creates bedfiles corresponding to 5 distinct regions:
# CDS, 3'UTR, 5'UTR, proximal intron (500bp), distal intron (500bp).
# Warning: this takes a long time!

# these two are really a minimum

from __future__ import print_function
from __future__ import division

# uncomment from this compatibility import list, as py3/py2 support progresses

from __future__  import absolute_import
from __future__  import unicode_literals
from future import standard_library
# from future.builtins import builtins
from future.builtins import utils
from future.utils import raise_with_traceback
from future.utils import iteritems

from argparse import ArgumentParser
import sys
import gffutils
import pybedtools
import pandas as pd
from collections import defaultdict
from . import annotation_functions as af
from tqdm import trange

def create_bedtools(features, keys, by_transcript=False):
    """Convert inclusive GFF features to half-open BED intervals.

    :param features: Iterable of one-based inclusive GFF features.
    :param keys: Source vocabulary from annotation_functions.get_keys.
    :param by_transcript: True to name intervals by transcript and preserve individual features; False to name by gene and merge exported overlapping intervals.
    :return: A BedTool with chromosome, converted start/end, selected ID, score 0, and strand. Multiple IDs produce separate intervals.
    """
    intervals = []
    key = 'transcript_id' if by_transcript else 'gene_id'
    progress = trange(len(features), desc='creating bedtools.')
    for feature in features:
        for i in range(len(feature.attributes[keys[key]])):
            interval = pybedtools.create_interval_from_list([
                feature.seqid, str(feature.start - 1), str(feature.end),
                feature.attributes[keys[key]][i], '0',
                feature.strand
            ])
            intervals.append(interval)
        progress.update(1)
    bedtool = pybedtools.BedTool(intervals)
    return bedtool


def merge_bedtool_by_gene(bedtool):
    """Merge BED intervals separately by chromosome, gene name, and strand.

    :param bedtool: BED6 BedTool; intervals within each chromosome and gene are merged strand-specifically.
    :return: A DataFrame with chrom/start/end/name/score/strand columns, sorted by chromosome and start. Requires the bedtools executable.
    """
    df = bedtool.to_dataframe()
    merged = pd.DataFrame(
        columns=['chrom','start','end','name','score','strand']
    )
    progress = trange(len(set(df['chrom'])))
    for chrom in set(df['chrom']):
        progress.set_description("Merging chromosome {}".format(chrom))
        dx = df[df['chrom']==chrom]
        dx.sort_values('start', inplace=True)
        dy = dx.groupby(['name']).apply(
            lambda x: pybedtools.BedTool.from_dataframe(x).merge(
                s=True, c='4,5,6', o='distinct,distinct,distinct'
            )
        )
        for d in dy:
            merged = pd.concat([merged, d.to_dataframe()], axis=0)
        progress.update(1)
    # Conversion from BedTool to DataFrame messes up the column order
    # Re-arrange here
    merged.columns = ['chrom','start','end','name','score','strand']
    merged.sort_values(by=['chrom','start'], inplace=True)
    return merged

def create_cds_region_bedfile(db, keys, cds_out, by_transcript=False):
    """Return CDS intervals and optionally write their selected-name BED export.

    :param db: gffutils.FeatureDB containing the requested feature types.
    :param keys: Source vocabulary from annotation_functions.get_keys.
    :param cds_out: Optional CDS BED output path, or None.
    :param by_transcript: True to name intervals by transcript and preserve individual features; False to name by gene and merge exported overlapping intervals.
    :return: The original unmerged CDS BedTool; a written export is merged by gene unless by_transcript is True.
    """
    cds_features = []

    for cds_feature in db.features_of_type(keys['cds']):
        cds_features.append(cds_feature)
    cds = create_bedtools(cds_features, keys, by_transcript)

    if by_transcript:
        cdsdf = cds.to_dataframe()
    else:
        cdsdf = merge_bedtool_by_gene(cds)

    if cds_out is not None:
        cdsdf.sort_values(by=['chrom','start','end'], inplace=True)
        cdsdf.to_csv(
            cds_out, sep=str("\t"), header=False,
            index=False
        )
    return cds

def create_exon_region_bedfile(db, keys, exon_out, by_transcript=False):
    """Write exon intervals as a BED export.

    :param db: gffutils.FeatureDB containing the requested feature types.
    :param keys: Source vocabulary from annotation_functions.get_keys.
    :param exon_out: Exon BED output path.
    :param by_transcript: True to name intervals by transcript and preserve individual features; False to name by gene and merge exported overlapping intervals.
    :return: None; writes a BED6 export, sorted when by_transcript is True.
    """
    exon_features = []

    for exon_feature in db.features_of_type(keys['exon']):
        exon_features.append(exon_feature)
    exons = create_bedtools(exon_features, keys, by_transcript)
    if not by_transcript:
        merge_bedtool_by_gene(exons).to_csv(
            exon_out, sep=str("\t"), header=False,
            index=False
        )
    else:
        exons = exons.sort()
        exons.to_dataframe().to_csv(
            exon_out, sep=str("\t"), header=False,
            index=False
        )


def create_gene_region_bedfile(db, keys, gene_out):
    """Write gene features as half-open BED intervals.

    :param db: gffutils.FeatureDB containing the requested feature types.
    :param keys: Source vocabulary from annotation_functions.get_keys.
    :param gene_out: Gene BED output path.
    :return: None; writes the unmerged gene intervals.
    """
    gene_features = []
    gene_hash = defaultdict(list)
    for gene_feature in db.features_of_type(keys['gene']):
        gene_features.append(gene_feature)
        #  TODO: remove this block of code when we know gene ids are unique
        for gene_id in gene_feature.attributes['gene_id']:
            if gene_id in gene_hash.keys():
                print("hashkey exists for gene {} at: {}".format(gene_id, gene_hash[gene_id]))
            gene_hash[gene_id].append(gene_feature.start)

    genes = create_bedtools(gene_features, keys)
    # TODO: double check this to make sure we don't need to merge by gene
    # (if more than 2 lines in the gtf file defines 1 gene, need to merge)
    genes.to_dataframe().to_csv(
        gene_out, sep=str("\t"), header=False,
        index=False
    )

    # merge_bedtool_by_gene(genes).to_csv(
    #     gene_out, sep=str("\t"), header=False,
    #     index=False
    # )

def create_utr_region_bedfiles(db, keys, cds_dict, utr3_out, utr5_out, by_transcript=False):
    """Return and optionally export five-prime and three-prime UTR intervals.

    :param db: gffutils.FeatureDB containing the requested feature types.
    :param keys: Source vocabulary from annotation_functions.get_keys.
    :param cds_dict: Mapping of transcript IDs to inclusive low/hi CDS bounds, used to classify generic UTR features.
    :param utr3_out: Optional three-prime UTR output path, or None.
    :param utr5_out: Optional five-prime UTR output path, or None.
    :param by_transcript: True to name intervals by transcript and preserve individual features; False to name by gene and merge exported overlapping intervals.
    :return: A tuple (five_prime, three_prime) of unmerged BedTools. Supports generic UTR or explicit utr5/utr3 feature types. Returns 1 with a diagnostic if neither vocabulary is available.
    """
    five_prime_utr_features = []
    three_prime_utr_features = []
    # Use CDS dict to determine for each UTR, whether it's 3/5'
    if keys['utr'] is not None:
        for utr_feature in db.features_of_type(keys['utr']):
            classified_utr = af.classify_utr(utr_feature, cds_dict)
            if classified_utr == '5utr':
                five_prime_utr_features.append(utr_feature)
            elif classified_utr == '3utr':
                three_prime_utr_features.append(utr_feature)
    elif keys['utr3'] is not None and keys['utr5'] is not None:
        for utr_feature in db.features_of_type(keys['utr5']):
            five_prime_utr_features.append(utr_feature)
        for utr_feature in db.features_of_type(keys['utr3']):
            three_prime_utr_features.append(utr_feature)
    else:
        print("Having trouble guessing the UTR key dictionary")
        return 1

    # print('key: {}, key: {}'.format(keys['utr3'], keys['utr5']))
    # Transform list of Intervals into BedTool
    utr5 = create_bedtools(five_prime_utr_features, keys, by_transcript)
    utr3 = create_bedtools(three_prime_utr_features, keys, by_transcript)
    # Merge any overlapping introns by transcript and save

    if by_transcript:
        utr5df = utr5.to_dataframe()
    else:
        utr5df = merge_bedtool_by_gene(utr5)

    if utr5_out is not None:
        utr5df.to_csv(
            utr5_out, sep=str("\t"), header=False, index=False
        )


    if by_transcript:
        utr3df = utr3.to_dataframe()
    else:
        utr3df = merge_bedtool_by_gene(utr3)

    if utr3_out is not None:
        utr3df.to_csv(
            utr3_out, sep=str("\t"), header=False, index=False
        )
    return utr5, utr3

def create_intron_region_bedfiles(db, exons_dict, transcripts_dict, keys,
            proxintron_out, distintron_out, allintron_out, by_transcript):
    """Return and optionally export proximal, distal, and complete intron intervals.

    :param db: gffutils.FeatureDB containing the requested feature types.
    :param exons_dict: Mapping of transcript IDs to inclusive exon coordinate dictionaries.
    :param transcripts_dict: Mapping of transcript IDs to inclusive start/end dictionaries.
    :param keys: Source vocabulary from annotation_functions.get_keys.
    :param proxintron_out: Optional proximal intron output path, or None.
    :param distintron_out: Optional distal intron output path, or None.
    :param allintron_out: Optional complete intron output path, or None.
    :param by_transcript: True to name intervals by transcript and preserve individual features; False to name by gene and merge exported overlapping intervals.
    :return: A tuple (proximal, distal, all_introns) of unmerged BedTools. Proximal regions lie within 500 bases of an exon; introns of length at most 1000 are entirely proximal.
    """
    allintrons = []
    proxintrons = []
    distintrons = []
    distance = 500

    # if we're finding all introns on a per-transcript bases, every name should contain the origin transcript id.
    # if not, we're going to collapse (merge) transcripts into nonoverlapping genic regions, so every name
    # will contain the gene_id.
    if by_transcript:
        name = keys['transcript_id']
    else:
        name = keys['gene_id']
    progress = trange(len(af.chromosome_set(db)))

    for chrom in af.chromosome_set(db):
        progress.set_description("Calculating introns for chromosome {}".format(chrom))

        for element in db.region(seqid=chrom):

            # make sure we JUST get the transcripts, gene features don't have associated exons
            if element.featuretype == 'transcript':
                # transcript_id is iterable here, one element may associate w/ multiple txids?
                for transcript_id in element.attributes[keys['transcript_id']]:
                    exons = exons_dict[transcript_id]
                    transcript = transcripts_dict[transcript_id]
                    # Given a transcript ID and its exons, find introns
                    introns = af.find_introns(transcript, exons)
                    # gene_id is iterable here, maybe one transcript may have multiple geneids?
                    for intron in introns:
                        for i in range(len(element.attributes[name])):
                            # Transform Feature into BedTool.Interval
                            # Important! Interval name is by gene, so outputs will
                            # merge all transcript introns by their gene id!
                            intron_interval = pybedtools.create_interval_from_list(
                                [element.chrom, str(int(intron['start']) - 1),
                                 intron['end'],  # turn 1 based into 0 based
                                 element.attributes[name][i],
                                 str(element.score), element.strand]
                            )
                            allintrons.append(intron_interval)
                            proxdist_dict = af.get_proxdist_from_intron(
                                intron_interval, distance=distance
                            )
                            for prox in proxdist_dict['prox']:
                                prox.name = element.attributes[name][i]
                                proxintrons.append(prox)
                            for dist in proxdist_dict['dist']:
                                dist.name = element.attributes[name][i]
                                distintrons.append(dist)
        progress.update(1)
    # Transform list of Intervals into BedTool
    proxintrons = pybedtools.BedTool(proxintrons)
    distintrons = pybedtools.BedTool(distintrons)
    allintrons = pybedtools.BedTool(allintrons)
    # Merge any overlapping introns by transcript and save

    if by_transcript:
        proxintronsdf = proxintrons.to_dataframe()
    else:
        proxintronsdf = merge_bedtool_by_gene(proxintrons)

    if proxintron_out is not None:
        proxintronsdf.to_csv(
            proxintron_out, sep=str("\t"), header=False, index=False
        )


    if by_transcript:
        distintronsdf = distintrons.to_dataframe()
    else:
        distintronsdf = merge_bedtool_by_gene(distintrons)

    if distintron_out is not None:
        distintronsdf.to_csv(
            distintron_out, sep=str("\t"), header=False, index=False
        )

    if by_transcript:
        allintronsdf = distintrons.to_dataframe()
    else:
        allintronsdf = merge_bedtool_by_gene(allintrons)

    if allintron_out is not None:
        allintronsdf.to_csv(
            allintron_out, sep=str("\t"), header=False, index=False
        )
    return proxintrons, distintrons, allintrons

def create_region_bedfiles(
        db_file, species, cds_out, utr3_out, utr5_out,
        proxintron_out, distintron_out, allintron_out, gene_out, exon_out, by_transcript
):
    ### get the parsed regions ###

    """Write selected BED exports from a feature database.

    :param db_file: Path to a gffutils database.
    :param species: Source vocabulary accepted by get_keys.
    :param cds_out: Optional CDS BED destination.
    :param utr3_out: Optional three-prime UTR BED destination.
    :param utr5_out: Optional five-prime UTR BED destination.
    :param proxintron_out: Optional proximal intron BED destination.
    :param distintron_out: Optional distal intron BED destination.
    :param allintron_out: Optional complete intron BED destination.
    :param gene_out: Optional gene BED destination.
    :param exon_out: Optional exon BED destination.
    :param by_transcript: Preserve transcript names instead of gene merging.
    :return: None; only non-None destinations select exports. Requires bedtools.
    """
    keys = af.get_keys(species)
    db = gffutils.FeatureDB(db_file)
    featuretypes = [x.lower() for x in list(db.featuretypes())]

    # Creates the CDS regions
    if cds_out is not None:
        create_cds_region_bedfile(db, keys, cds_out, by_transcript)

    # Creates 3 and 5' UTR regions
    if utr3_out is not None or utr5_out is not None:
        cds_dict = af.get_all_cds_dict(db, keys['cds'])
        create_utr_region_bedfiles(db, keys, cds_dict, utr3_out, utr5_out, by_transcript)

    # Creates the prox and dist intron regions
    if proxintron_out is not None or distintron_out is not None or allintron_out is not None:
        exons_dict = af.get_all_exons_dict(db, keys['exon'], keys['transcript_id'])
        transcripts_dict = af.get_all_transcripts_dict(db, keys['exon'], keys[
            'transcript_id'])
        create_intron_region_bedfiles(
            db, exons_dict, transcripts_dict, keys,
            proxintron_out, distintron_out, allintron_out, by_transcript
        )
    # Creates the gene regions
    if gene_out is not None:
        create_gene_region_bedfile(db, keys, gene_out)

    # Creates the exon regions
    if exon_out is not None:
        create_exon_region_bedfile(db, keys, exon_out, by_transcript)

def main():
    # Setup argument parser
    """Export requested BED regions from command-line db_file and species options.

    Reads sys.argv. Optional cds/utr3/utr5/proxintron/distintron/allintron/gene/exon
    output paths select exports; by-transcript retains transcript names.
    Returns None after writing requested files. Requires bedtools.
    """
    parser = ArgumentParser()

    parser.add_argument(
        "--db_file",
        dest="db_file",
        help="gtfdb file",
        required=True
    )
    parser.add_argument(
        "--species",
        dest="species",
        help="sets the species gtf nomenclature",
        required=True
    )
    parser.add_argument(
        "--cds_out",
        dest='cds_out',
        help="output of cds bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--utr5_out",
        dest='utr5_out',
        help="output of utr5 bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--utr3_out",
        dest='utr3_out',
        help="output of utr3 bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--proxintron_out",
        dest='proxintron_out',
        help="output of proxintron bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--distintron_out",
        dest='distintron_out',
        help="output of distintron bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--allintron_out",
        dest='allintron_out',
        help="output of all introns bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--gene_out",
        "--genes_out",
        dest='gene_out',
        help="output of genes bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--exon_out",
        "--exons_out",
        dest='exon_out',
        help="output of exons bed",
        required=False,
        default=None
    )
    parser.add_argument(
        "--by-transcript",
        "--by_transcript",
        dest='by_transcript',
        help="output CDS/UTR3/UTR5 by transcript (Don't do gene merge)",
        required=False,
        default=False,
        action='store_true'
    )
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    db_file = args.db_file
    cds_out = args.cds_out
    utr3_out = args.utr3_out
    utr5_out = args.utr5_out
    proxintron_out = args.proxintron_out
    distintron_out = args.distintron_out
    allintron_out = args.allintron_out
    gene_out = args.gene_out
    exon_out = args.exon_out
    species = args.species
    by_transcript = args.by_transcript

    create_region_bedfiles(db_file, species, cds_out, utr3_out, utr5_out,
                           proxintron_out, distintron_out, allintron_out,
                           gene_out, exon_out, by_transcript)
if __name__ == "__main__":
    main()