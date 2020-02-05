import pandas as pd
import numpy as np
import gffutils
import sys
import pybedtools
from collections import defaultdict
from future.utils import iteritems

MAXVAL = 1000000000
MINVAL = 0


def get_keys(species):
    """
    Describes the language for each kind of gtf file 
    (usually different between sources).
    Defaults to gencode standard for hg19/mm9/mm10/hg38

    :param species: string
        either one of: 'hg19_gencode','mm10','ce11','mm9','hg38'
    :return:
    """
    print("species is {}".format(species))
    if species == 'hg19_ensembl':
        cds_key = 'CDS'
        utr3_key = 'three_prime_utr'
        utr5_key = 'five_prime_utr'
        utr_key = None
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    elif species == 'ce10':
        cds_key = 'CDS'
        utr3_key = None # 'three_prime_UTR'
        utr5_key = None # 'five_prime_UTR'
        utr_key = 'UTR'
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    elif species == 'ce11':
        cds_key = 'CDS'
        utr3_key = 'three_prime_utr'
        utr5_key = 'five_prime_utr'
        utr_key = None
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    else:
        cds_key = 'CDS'
        utr3_key = None  #
        utr5_key = None  # in human/mice, this key doesn't exist
        utr_key = 'UTR'
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_type'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_type'

    keys = {
        'cds': cds_key,
        'utr3': utr3_key,
        'utr5': utr5_key,
        'utr': utr_key,
        'gene': gene_key,
        'gene_name': gene_name_key,
        'transcript': transcript_key,
        'transcript_id': transcript_id_key,
        'transcript_type': type_key,
        'exon': exon_key,
        'gene_id': gene_id_key,
        'gene_type': gene_type_key
    }
    return keys


def get_gene_to_transcript_dict(db, gene_id_key, transcript_id_key):
    """
    Returns a gene: transcript dictionary. 
    
    gene:[transcript1, transcript2]
    
    :param db: gffutils.FeatureDB
    :param gene_id_key: string
    :param transcript_id_key: string
    :return: 
    """

    genes_dict = defaultdict(list)
    for gene_feature in db.features_of_type('transcript'):
        for gene_id in gene_feature.attributes[gene_id_key]:
            for transcript_id in gene_feature.attributes[transcript_id_key]:
                genes_dict[gene_id].append(transcript_id)
    return genes_dict


def get_longest_transcripts(genes_dict, transcripts_dict):
    """
    Returns a dictionary of genes : longest_transcript
    """
    longest_genes_dict = defaultdict(dict)
    for gene, transcripts in iteritems(genes_dict):
        max_transcript_len = -1
        max_transcript = ""
        for transcript in transcripts:
            transcript_len = transcripts_dict[transcript]['end'] - transcripts_dict[transcript]['start']
            if transcript_len > max_transcript_len:
                max_transcript_len = transcript_len
                max_transcript = transcript
        longest_genes_dict[gene] = max_transcript
    return longest_genes_dict


def most_upstream_downstream_positions(genes_dict, transcripts_dict):
    """
    Assumes the most downstream start < most downstream end (+, flipped for -) foreach transcript.
    Returns a dictionary of genes : 
    """
    d = defaultdict(dict)
    for gene, transcripts in iteritems(genes_dict):
        min_transcript_pos = 1000000000    # as long as we don't have any chromosomes larger than 1 billion
        max_transcript_pos = -1
        for transcript in transcripts:
            if transcripts_dict[transcript]['end'] > max_transcript_pos:
                max_transcript_pos = transcripts_dict[transcript]['end']
            if transcripts_dict[transcript]['start'] < min_transcript_pos:
                min_transcript_pos = transcripts_dict[transcript]['start']

        d[gene] = {'start': min_transcript_pos, 'end': max_transcript_pos}
    return d

def get_premrna_lengths(gene_id, upstream_downstream):
    start = upstream_downstream[gene_id]['start']
    end = upstream_downstream[gene_id]['end']
    return end - start + 1    # (coords are 1-based inclusive still)


def get_mrna_lengths(gene_id, exons_dict, genes_dict):
    exons_list = []
    total_mrna_length = 0
    for transcript in genes_dict[gene_id]:
        for exon in exons_dict[transcript]:
            exons_list.append(pybedtools.create_interval_from_list(
                [exon['chrom'],
                 str(exon['start'] - 1),  # since we're converting to bedtool, use 0-based
                 str(exon['end']),
                 transcript,
                 '0',
                 exon['strand']]
            ))
    exons_list = pybedtools.BedTool(exons_list)
    exons_list = exons_list.sort().merge()  # not strand specific
    for exon in exons_list:
        total_mrna_length += (exon.end - exon.start)
    return total_mrna_length


def get_all_cds_dict(db, cds_key):
    """
    For every cds-annotated transcript id (ENST), return a
    dictionary containing the lowest and highest
    cds start and end vals for that transcript.

    :return chr19_cds_dict : defaultdict{transcript:{'start':START, 'end':END}}
    """
    # chr19_cds_dict = defaultdict(lambda: {'low': MAXVAL, 'hi': MINVAL})
    cds_dict = defaultdict(dict)
    for cds_feature in db.features_of_type(cds_key):
        for transcript_id in cds_feature.attributes['transcript_id']:
            # if cds_feature.start <= chr19_cds_dict[transcript_id]['low']:
            if cds_feature.start <= cds_dict[transcript_id].get("low", MAXVAL):
                cds_dict[transcript_id]['low'] = cds_feature.start
            # if cds_feature.end >= chr19_cds_dict[transcript_id]['hi']:
            if cds_feature.end >= cds_dict[transcript_id].get("hi", MINVAL):
                cds_dict[transcript_id]['hi'] = cds_feature.end
    return cds_dict


def get_all_exons_dict(db, exon_key, transcript_id_key):
    """
    Returns dictionary of exons as transcript_id:{
        [
            {'start':START, 'end':END},
            {'start':START, 'end':END},
            ...
        ]
    }.
  
    :param db: gffutils.FeatureDB
    :param transcript_id_key: string
    :param exon_key: string
    :return: 
    """

    exons_dict = defaultdict(list)
    for exon_feature in db.features_of_type(exon_key):
        for transcript_id in exon_feature.attributes[transcript_id_key]:
            exons_dict[transcript_id].append(
                {
                    'chrom': exon_feature.seqid,
                    'start': exon_feature.start,
                    'end': exon_feature.end,
                    'strand': exon_feature.strand,
                }
            )
    return exons_dict


def get_all_transcripts_dict(db, transcript_key, transcript_id_key):
    """
    Returns dictionary of transcript_id:{'start':START, 'end':END}.
  
    :param db: gffutils.FeatureDB
    :param transcript_key: string
    :param transcript_id_key: string
    :return: chr19_transcripts_dict: defaultdict(dict)
        hash of transcripts and their start/end coordinates
    """
    transcripts_dict = defaultdict(dict)
    for transcript_feature in db.features_of_type(transcript_key):
        for transcript_id in transcript_feature.attributes[transcript_id_key]:
            transcripts_dict[transcript_id] = {
                'start': transcript_feature.start,
                'end': transcript_feature.end
            }
    return transcripts_dict


def chromosome_set(db):
    """
    Returns the set of chromosomes that exist in a database.

    :return:
    """
    ret = db.execute("SELECT seqid FROM features").fetchall()
    all_chromosomes = [r['seqid'] for r in ret]
    return set(all_chromosomes)


def gene_id_to_name(db, gene_name_key):
    """
    Returns a dictionary containing a gene_id:name translation
    Note: may be different if the 'gene_id' or 'gene_name'
    chr19_keys are not in the source GTF file
    (taken from gscripts.region_helpers)

    :param db: gffutils.FeatureDB
        gffutils-created feature database (from 
        gffutils.create_db())
    :param gene_name_key: string
        Typically gene_name
    :return gene_name_dict : dict
        dict of {gene_id : gene_name}
    """

    genes = db.features_of_type('gene')
    gene_name_dict = {}

    for gene in genes:
        try:
            gene_id = gene.attributes['gene_id'][0] if type(
                gene.attributes['gene_id']
            ) == list else gene.attributes['gene_id']
            gene_name_dict[gene_id] = gene.attributes[gene_name_key][0]
        except KeyError:
            print(gene.attributes.keys())
            print("Warning. Key not found for {}".format(gene))
            sys.exit(1)
    return gene_name_dict


def classify_utr(utr_feature, cds_dict):
    """
    Given a feature classified as a UTR, return whether or not it is
    upstream (5') or downstream (3') based on CDS positions

    :param utr_feature: gffutils.Feature
        feature already classified as UTR (within a coding transcript but
        outside of CDS regions).
    :param cds_dict: dict
        dictionary containing cds positions for every transcript
    :return:
    """
    three_prime_utr = False
    five_prime_utr = False

    for transcript_id in utr_feature.attributes['transcript_id']:
        try:
            if utr_feature.strand == '+':
                if cds_dict[transcript_id]['low'] > utr_feature.end:
                    five_prime_utr = True
                if cds_dict[transcript_id]['hi'] < utr_feature.start + 1:
                    three_prime_utr = True
            elif utr_feature.strand == '-':
                if cds_dict[transcript_id]['low'] > utr_feature.end:
                    three_prime_utr = True
                if cds_dict[transcript_id]['hi'] < utr_feature.start + 1:
                    five_prime_utr = True
        except KeyError as e:
            return 'unclassified_utr' # no CDS found for this transcript so cannot positionally assign utr
    if five_prime_utr:
        return '5utr'
        # return 'five_prime_utr'
    elif three_prime_utr:
        return '3utr'
        # return 'three_prime_utr'
    else:
        return 'unclassified_utr'


def find_introns(transcript, exons):
    """
    Given transcript coordinates and a list of exons,
    return a list of introns (inverse coordinates of exons)
    This is 1-based inclusive, given 1-based inclusive exon coords

    :param transcript: dictionary
        {'start':int, 'end':int}
    :param exons: list[dict]
        [
            {'start':int,'end':int},
            {'start':int,'end':int},
            ...
        ]
    :return introns: list[dict]
        [
            {'start':int,'end':int},
            {'start':int,'end':int},
            ...
        ]
    """
    positions = [] # list of
    introns = []

    pos_append = positions.append
    intron_append = introns.append

    for exon in exons:
        pos_append(exon['start'] - 1)
        pos_append(exon['end'] + 1)
    positions = sorted(positions)
    try:
        # feature doesn't start with an intron
        if positions[0] < transcript['start']:
            positions.pop(0)
        else:
            positions.insert(transcript['start'])
        # feature doesn't end at an intron
        if positions[-1] > transcript['end']:
            positions.pop(-1)
        else:
            pos_append(transcript['end'])
        for i in range(0, len(positions) - 1, 2):
            intron_append({'start': positions[i], 'end': positions[i + 1]})
    except IndexError:
        print("No exons found: {}".format(transcript), positions)
        return []
    return introns


def get_proxdist_from_intron(interval, distance=500):
    """
    Given an interval, return a list of prox and dist intervals.
    Note: This uses pybedtools, which is 0 based. So unlike 1-based
    gffutils which is used in most of these functions, it returns
    zero-based exclusive bedtools intervals!

    :param interval: pybedtools.Interval
        intron interval
    :param distance: int
        maximum distance away from an exon that an intron position can be
        before calling it a distal intron region.
    :return proxdist: dict(list)
        dictionary containing 'prox' and 'dist' as chr19_keys, list of
        each proximal and distal intron region as values.
    """

    d = distance * 2 # total distance

    # if the region is less than 2x the distance, any point on
    # that region will be equal or less than an exon.
    if interval.end - interval.start <= d:
        prox = [
            interval.chrom,
            '{}'.format(interval.start),
            '{}'.format(interval.end),
            'proxintron{}'.format(distance),
            '{}'.format(interval.score),
            interval.strand
        ]
        prox = pybedtools.create_interval_from_list(prox)
        return {'prox': [prox], 'dist': []}
    else:
        prox_left = [
            interval.chrom,
            '{}'.format(interval.start),
            '{}'.format(interval.start + distance),
            'proxintron{}'.format(distance),
            '{}'.format(interval.score),
            interval.strand
        ]
        prox_left = pybedtools.create_interval_from_list(prox_left)
        prox_right = [
            interval.chrom,
            '{}'.format(interval.end - distance),
            '{}'.format(interval.end),
            'proxintron{}'.format(distance),
            '{}'.format(interval.score),
            interval.strand
        ]
        prox_right = pybedtools.create_interval_from_list(prox_right)
        dist = [
            interval.chrom,
            '{}'.format(interval.start + distance),
            '{}'.format(interval.end - distance),
            'distintron{}'.format(distance),
            '{}'.format(interval.score),
            interval.strand
        ]
        dist = pybedtools.create_interval_from_list(dist)
        return {'prox': [prox_left, prox_right], 'dist': [dist]}