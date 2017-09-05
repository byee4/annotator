#!/usr/bin/env python

# transitionning to python2/python3 support

# these two are really a minimum

from __future__ import print_function
from __future__ import division

# uncomment from this compatibility import list, as py3/py2 support progresses

# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
from future.utils import iteritems

# import pybedtools
from tqdm import trange
import gffutils
from collections import defaultdict
from collections import OrderedDict
import copy
import sys
import csv
import multiprocessing
import concurrent.futures
from itertools import chain

import cProfile

HASH_VAL = 1000000
MAXVAL = 1000000000
MINVAL = 0
FEATURES_DICT = {}

def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func

def create_definitions(db_file, chroms=[], species='hg19', append_chr=False, fuzzy=0):
    """
    This method creates definitions and populates dictionaries to be used
    for annotation. Previously initialized in Annotator().

    :param db_file:
    :param chroms:
    :param species:
    :param append_chr:
    :param fuzzy:
    :return:
    """

    cds_key, utr3_key, utr5_key, utr_key, \
    gene_name_key, transcript_id_key, type_key = get_keys(species)

    db = gffutils.FeatureDB(db_file)
    featuretypes = [x.lower() for x in list(db.featuretypes())]

    progress = trange(5, desc='Initializing/creating defs')
    progress.set_description("Building gene ids -> gene names dictionary")

    geneid_to_name_dict = gene_id_to_name(db, gene_name_key)
    progress.update(1)
    if len(chroms) != 0:  # if use specific chromosomes, otherwise hash all
        chromosomes = set(chroms)
    else:
        chromosomes = chromosome_set(db)

    exons_dict = {}
    transcripts_dict = {}
    cds_dict = {}

    ### If 'exon' features exist in db, catalog it.
    if 'exon' in featuretypes:
        progress.set_description("Adding all exon boundaries")
        exons_dict = get_all_exons_dict(db, transcript_id_key)
    progress.update(1)
    ### If 'transcript' features exist in db, catalog it.
    if 'transcript' in featuretypes:
        progress.set_description("Adding all transcript boundaries")
        transcripts_dict = get_all_transcripts_dict(db, transcript_id_key)
    progress.update(1)
    ### If 'cds' features exist in db, catalog it.
    if 'cds' in featuretypes:
        progress.set_description("Adding all CDS boundaries")
        cds_dict = get_all_cds_dict(db, cds_key)
    progress.update(1)
    progress.set_description("Adding all features in DB to dictionary")
    features_dict, num_features = hash_features(
        db, chromosomes, append_chr, fuzzy,
        exons_dict, transcripts_dict,
        transcript_id_key, featuretypes
    )
    progress.update(1)
    return geneid_to_name_dict, exons_dict, transcripts_dict, cds_dict, features_dict, \
           cds_key, utr3_key, utr5_key, utr_key, gene_name_key, transcript_id_key, type_key

def get_keys(species):
    """
    Describes the language for each kind of gtf file (usually different between sources).

    :param species:
    :return:
    """
    if species == 'ce11':
        cds_key = 'CDS'
        utr3_key = 'three_prime_UTR'
        utr5_key = 'five_prime_UTR'
        utr_key = None
        gene_name_key = 'gene_id'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
    else:
        cds_key = 'CDS'
        utr3_key = None  #
        utr5_key = None  # in human/mice, this key doesn't really exist
        utr_key = 'UTR'
        gene_name_key = 'gene_name'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_type'

    return cds_key, utr3_key, utr5_key, utr_key, gene_name_key, transcript_id_key, type_key

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
    keys are not in the source GTF file
    (taken from gscripts.region_helpers)

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

def hash_features(db, chromosomes, append_chr, fuzzy, exons_dict, transcripts_dict, transcript_id_key, featuretypes):
    """
    hashes features by position.
    :return features_dict : collections.defaultdict()
        dictionary of features{[chrom, pos/HASH_VAL, strand] : feature_list}
    """
    num_features = 0

    features_dict = defaultdict(list)

    progress = trange(
        len(chromosomes),
        leave=False,
        desc='Building feature index'
    )

    for chrom in chromosomes:
        fixed_chrom = 'chr' + chrom if append_chr else chrom

        progress.set_description("Adding {}...".format(chrom))

        for element in db.region(seqid=chrom):
            introns = []

            if element.featuretype == 'gene':
                if element.strand == '+':
                    element.start = element.start - fuzzy
                elif element.strand == '-':
                    element.end = element.end + fuzzy
            if element.featuretype == 'gene':
                if element.strand == '+':
                    element.end = element.end + fuzzy
                elif element.strand == '-':
                    element.start = element.start - fuzzy
            if element.featuretype == 'transcript' and 'intron' not in featuretypes: # if we see at least one intron, we assume introns are annotated!
                for transcript_id in element.attributes[transcript_id_key]:
                    exons = exons_dict[transcript_id]
                    transcript = transcripts_dict[transcript_id]
                    introns = find_introns(transcript, exons)

            start = int(element.start / HASH_VAL)
            end = int(element.end / HASH_VAL)

            for i in range(start, end + 1):
                element.chrom = fixed_chrom
                append = features_dict[fixed_chrom, i, element.strand].append
                append(element)

                for intron in introns:
                    intron_feature = copy.deepcopy(element)
                    intron_feature.start = intron['start']
                    intron_feature.end = intron['end']
                    intron_feature.featuretype = 'intron'
                    append(
                        intron_feature
                    )
            num_features += 1
        progress.update(1)

    print("NUM FEATURES: {}".format(num_features))
    return features_dict, num_features

def get_all_exons_dict(db, transcript_id_key):
    """
    Returns dictionary of exons as transcript_id:{
        [
            {'start':START, 'end':END},
            {'start':START, 'end':END},
            ...
        ]
    }.

    :return:
    """
    exons_dict = defaultdict(list)
    for exon_feature in db.features_of_type('exon'):
        for transcript_id in exon_feature.attributes[transcript_id_key]:
            exons_dict[transcript_id].append(
                {
                    'start': exon_feature.start,
                    'end': exon_feature.end
                }
            )
    return exons_dict

def get_all_transcripts_dict(db, transcript_id_key):
    """
    Returns dictionary of transcript_id:{'start':START, 'end':END}.

    :return transcripts_dict: defaultdict(dict)
        hash of transcripts and their start/end coordinates
    """
    transcripts_dict = defaultdict(dict)
    for transcript_feature in db.features_of_type('transcript'):
        for transcript_id in transcript_feature.attributes[transcript_id_key]:
            transcripts_dict[transcript_id] = {
                'start': transcript_feature.start,
                'end': transcript_feature.end
            }
    return transcripts_dict


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
    positions = []
    introns = []
    pos_append = positions.append
    intron_append = introns.append

    for exon in exons:
        pos_append(exon['start'] - 1)
        pos_append(exon['end'] + 1)
    positions = sorted(positions)

    if positions[0] < transcript[
        'start']:  # there is no intron at the start of the feature
        positions.pop(0)
    else:
        positions.insert(transcript['start'])
    if positions[-1] > transcript['end']:
        positions.pop(-1)
    else:
        pos_append(transcript['end'])
    for i in range(0, len(positions) - 1, 2):
        intron_append({'start': positions[i], 'end': positions[i + 1]})
    return introns

def get_all_cds_dict(db, cds_key):
    """
    For every cds-annotated transcript id (ENST), return a
    dictionary containing the lowest and highest
    cds start and end vals for that transcript.

    :return cds_dict : defaultdict{transcript:{'start':START, 'end':END}}
    """
    # cds_dict = defaultdict(lambda: {'low': MAXVAL, 'hi': MINVAL})
    cds_dict = defaultdict(dict)
    for cds_feature in db.features_of_type(cds_key):
        for transcript_id in cds_feature.attributes['transcript_id']:
            # if cds_feature.start <= cds_dict[transcript_id]['low']:
            if cds_feature.start <= cds_dict[transcript_id].get("low",MAXVAL):
                cds_dict[transcript_id]['low'] = cds_feature.start
            # if cds_feature.end >= cds_dict[transcript_id]['hi']:
            if cds_feature.end >= cds_dict[transcript_id].get("hi",MINVAL):
                cds_dict[transcript_id]['hi'] = cds_feature.end
    return cds_dict

def classify_utr(utr_feature, cds_dict):
    """
    Given a list of features, return a dictionary of
    gene_id : {start: cds_start, end: cds_end} sites

    :param features : list[gffutils.Feature]
    :return cds_ends : dict
    """
    three_prime_utr = False
    five_prime_utr = False

    for transcript_id in utr_feature.attributes['transcript_id']:
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

    if five_prime_utr:
        return '5utr'
        # return 'five_prime_utr'
    elif three_prime_utr:
        return '3utr'
        # return 'three_prime_utr'
    else:
        return 'unclassified_utr'

def get_all_overlapping_features_from_query_stranded(
        chrom, qstart, qend, strand, features_dict):
    """
    Given a query location (chr, start, end), return all features that
    overlap by at least one base. Functions similarly to gffutils db.region(),
    but uses the pre-hashed self.features_dict to greatly speed things up.

    :param chrom : string
    :param qstart : int
    :param qend : int
    :param strand : string
    :return features : list
        list of gffutils.Feature objects.
    """
    features = []
    append = features.append
    start_key = int(qstart / HASH_VAL)
    end_key = int(qend / HASH_VAL)
    qstart = qstart + 1  # change 0-based bed to 1-based gff
    for i in range(start_key, end_key + 1):
        for feature in features_dict[chrom, i, strand]:
            if qstart <= feature.start and qend >= feature.end:  # query completely contains feature
                feature.attributes['overlap'] = 'query_contains_feature'
                append(feature)
            elif qstart >= feature.start and qend <= feature.end:  # feature completely contains query
                feature.attributes['overlap'] = 'feature_contains_query'
                append(feature)
            elif qstart <= feature.start and qend > feature.start:  # feature partially overlaps (qstart < fstart < qend)
                feature.attributes[
                    'overlap'] = 'partial_by_{}_bases'.format(
                    qend - (feature.start - 1))
                append(feature)
            elif qstart <= feature.end and qend >= feature.end:  # feature partially overlaps (qstart < fend < qend)
                feature.attributes[
                    'overlap'] = 'partial_by_{}_bases'.format(
                    feature.end - qstart)
                append(feature)
    # for f in features:
    #     f.start = f.start - 1 # fix to make gff-coordinates 0-based
    return features

def get_all_overlapping_features_from_query_unstranded(
        chrom, qstart, qend, strand, features_dict):
    """
    Returns overlapping features from a query.

    :param chrom:
    :param qstart:
    :param qend:
    :param strand:
    :return:
    """
    if strand is None:  # If no strand given, return all features from both strands
        positive_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, '+', features_dict
        )
        negative_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, '-', features_dict
        )
        return positive_features + negative_features
    else:  # return features on the same strand, only return features on opposite strand if it can't find anything otherwise.
        priority_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, strand, features_dict
        )
        if len(
                priority_features) == 0:  # found nothing on the proper strand
            if strand == '+':
                return get_all_overlapping_features_from_query_stranded(
                    chrom, qstart, qend, '-', features_dict
                )
            elif strand == '-':
                return get_all_overlapping_features_from_query_stranded(
                    chrom, qstart, qend, '+', features_dict
                )
            else:
                print("Strand not correct: {}".format(strand))
                return []
        else:
            return priority_features

def get_all_overlapping_features_from_query(
        chrom, qstart, qend, strand, features_dict, stranded=True):
    """
    Returns overlapping features from a query.
    If stranded is True, then return only the features that overlap on the same strand.
    Otherwise, features returned need not to be with respect to strand (although if
    strand is given, it will return same-stranded features first before attempting to
    find those on the opposite strand).

    :param chrom:
    :param qstart:
    :param qend:
    :param strand:
    :param stranded:
    :return:
    """
    if stranded:
        return get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, strand, features_dict
        )
    else:
        return get_all_overlapping_features_from_query_unstranded(
            chrom, qstart, qend, strand, features_dict
            # TODO: implement the priority better. For now, this will return all features for both strands
        )

def parse_annotation_string(features_string):
    """
    Splits a feature string into a list of feature strings

    :param features_string:
    :return:
    """
    features = features_string.split('|')
    return features

def is_protein_coding(transcript_type):
    """
    if defined protein coding, return True else False

    :param transcript_type:
    :return:
    """
    if transcript_type == 'protein_coding':
        return True

    return False

def return_highest_priority_feature(formatted_features, priority):
    """
    Given a list of formatted features, return the feature with the highest
    priority. If there are ties, return just the first one.

    :param formatted_features: list
        list of feature_strings
    :param priority: list
        list containing feature types in order of priority.
    :return:
    """
    # Build dict
    combined_dict = defaultdict(list)
    for feature_string in formatted_features:
        transcript, start, end, strand, feature_type, gene_id, \
        gene_name, transcript_type_list, overlap = feature_string.split(':')
        transcript_type_list = transcript_type_list.split(',')
        for transcript_type in transcript_type_list:
            if is_protein_coding(
                    transcript_type):  # simplify all the types at first
                combined_dict['protein_coding', feature_type].append(
                    feature_string)
            else:
                combined_dict['non_coding',  feature_type].append(
                    feature_string)
    # return the highest one
    combined_dict = OrderedDict(
        combined_dict
    )  # turn into ordered dict, is that ok?
    # print("combined dict: {}".format(combined_dict))
    try:
        # print("PRIORITY: ", priority)
        for key, _ in iteritems(combined_dict): # append nonexisting regions to the end of the priority list
            if list(key) not in priority:
                priority.append(list(key))
        # print("PRIORITY: ",priority)
        combined_dict = sorted(  # sort based on priority list
            iteritems(combined_dict),
            key=lambda x: priority.index([x[0][0], x[0][1]])
        )
        # print("PRIORITY: {}".format(priority))
        # print("combined dict (sorted): {}".format(combined_dict))
    except ValueError as e:
        print(e)
        print(combined_dict)
    return combined_dict[0]

def prioritize_transcript(unique_transcript_features, transcript_priority):
    """
    Given a dictionary of features, return a dictionary containing
    a singular transcript for each unique gene in the feature set.

    :param unique_transcript_features: defaultdict[list]
        dictionary where keys = transcript ids and values = list of
        features associated with each transcript that overlaps with a query
        (ie. {'ENSTXYZ':['CDS', 'exon', 'transcript']}
    :param transcript_priority: list[tuple]
        list of tuples (featuretype, transcripttype)
    :return unique_genes: defaultdict[list]
    """
    unique_transcripts = defaultdict(list)
    unique_genes = defaultdict(list)


    if len(unique_transcript_features.keys()) == 0:
        return 'intergenic'

    ### PRIORITIZE TRANSCRIPT ###
    for transcript in unique_transcript_features.keys():  # For each unique transcript
        top_transcript = return_highest_priority_feature(
            unique_transcript_features[transcript],
            transcript_priority
        )[1][0]  # [0] contains the dictionary key
        # unique_transcripts[transcript].append(
        #     top_transcript
        # )
        # add gene key
        gene_list = top_transcript.split(':')[5].split(',')
        for gene in gene_list:
            unique_genes[gene].append(top_transcript)
    return unique_genes

def prioritize_genes(unique_genes, gene_priority):
    final_transcripts = []
    final_transcripts_append = final_transcripts.append
    for gene, transcripts in iteritems(unique_genes):
        for transcript in transcripts:
            final_transcripts_append(transcript)

    feature_type, final_genes = return_highest_priority_feature(
        final_transcripts, gene_priority
    )

    if feature_type[0] == 'non_coding':  # TODO: fix. hacky
        ff = []
        for f in sorted(final_genes):
            f = f.replace(
                'exon', 'noncoding_exon'
            ).replace(
                'intron', 'noncoding_intron'
            )
            ff.append(f)
        return ff
    else:
        return sorted(final_genes)

def prioritize_transcript_then_gene(
        formatted_features, transcript_priority, gene_priority, parent='gene'):
    """
    Given a list of features, group them first by transcript
    and return the highest priority feature for each transcript.
    Then group each transcript into associated genes and
    return the highest priority feature.

    :param formatted_features: list[string]
        list of features as string representation
    :param transcript_priority: list
        list of tuples (featuretype, transcripttype)
    :param gene_priority: list
        list of tuples (featuretype, transcripttype
    :param parent: string
        the highest/most general annotation in the schema (usually gene).
    :return top_gene: string
        string representation of the highest priority feature
        after sorting by transcript priority, then by
        gene priority
    """
    unique_transcripts = defaultdict(list)

    for feature_string in formatted_features:
        if feature_string.split(':')[4] != parent:
            transcript = feature_string.split(':')[0]
            unique_transcripts[transcript].append(
                feature_string
            )

    unique_genes = prioritize_transcript(unique_transcripts, transcript_priority)
    top_gene = prioritize_genes(unique_genes, gene_priority)
    return top_gene

def split_string(priority_string, delimiter=':'):
    """
    Parses/splits a string with an expected format and returns
    three important fields: region, gene_id, gene_name. Returns
    'intergenic' for all three fields if intergenic

    :param priority_string: string
        colon-delimited string containing region geneid and
        genename in fields 4-6 (0-based).

    :return region, gene_id, gene_name:
    """
    if priority_string == 'intergenic':
        return 'intergenic', 'intergenic', 'intergenic'
    else:
        fields = priority_string.split(delimiter)
        return fields[4], fields[5], fields[6]

def annotate(
        chrom, start, end, name, score, strand,
        stranded, transcript_priority, gene_priority,
        features_dict, cds_dict, transcript_id_key, type_key):

    """
    Given an interval, annotates using the priority

    :param interval:
    :return:
    """
    append_count = 0

    overlapping_features = get_all_overlapping_features_from_query(
        chrom, start, end, strand, features_dict, stranded
    )

    #   If we find fno overlapping features, return 'intergenic' (no feature seen)
    if len(overlapping_features) == 0: # gene, name, region, type, to_append
        return chrom, start, end, name, score, strand, \
               'intergenic' , 'intergenic', \
               'intergenic', 'intergenic'
    to_append = ''  # full list of genes overlapping features
    transcript = defaultdict(list)
    for feature in overlapping_features:  # for each overlapping feature
        #   for M10 annotations, I see gene featuretypes without transcript_ids.
        #   That's fine. We can re-construct genes from the transcript features, so ignore
        #   anything that doesn't have a transcript_id in it.
        if transcript_id_key in feature.attributes.keys():
            for transcript_id in feature.attributes[
                transcript_id_key
            ]:  # multiple genes can be associated with one feature
                transcript[transcript_id].append(
                    feature)  # append features to their respective genes

    for transcript, features in iteritems(transcript):
        for feature in features:
            append_count += 1
            if feature.featuretype == 'UTR':  # if we only have UTR, specify what kind of UTR (3' or 5').
                feature.featuretype = classify_utr(feature, cds_dict)
            to_append += "{}:{}:{}:{}:{}:".format(
                transcript,
                feature.start - 1,  # report 0 based
                feature.end,
                feature.strand,
                feature.featuretype,
            )
            if 'gene_id' in feature.attributes.keys():
                for t in feature.attributes['gene_id']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            else:
                to_append = to_append + '-:'
            if 'gene_name' in feature.attributes.keys():
                for t in feature.attributes['gene_name']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            elif 'transcript_id' in feature.attributes.keys():
                for t in feature.attributes['transcript_id']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            else:
                to_append = to_append + '-:'
            if type_key in feature.attributes.keys():
                for t in feature.attributes[type_key]:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            else:
                to_append = to_append + '-:'
            if 'overlap' in feature.attributes.keys():
                for t in feature.attributes['overlap']:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + '|'
            else:
                to_append = to_append + '-:'

    to_append = to_append[:-1]

    if append_count == 1:  # if just one region, save computation by returning it.
        region, gene, rname = split_string(to_append)
    else:
        priority = prioritize_transcript_then_gene(
            parse_annotation_string(to_append),
            transcript_priority,
            gene_priority
        )
        region, gene, rname = split_string(priority[0])
        for highest_priority_gene in priority[1:]:
            current_region, current_gene, current_rname = split_string(highest_priority_gene)
            assert region == current_region  # there should really just be one highest priority region.
            ### if the highest priority region contains multiple genes, return all
            if gene != current_gene:
                gene = gene + ',' + current_gene
            if rname != current_rname:
                rname = rname + ',' + current_rname
    # print(feature.start, feature.end, feature.strand, to_append)
    return chrom, start, end, name, score, strand, gene , rname, region, to_append

def split_bed_interval(bed_string):
    """
    Takes a bed file as a line (tab delimited) and returns 6 fields
    correctly typed.
    :param bed_string:
    :return:
    """
    bed = bed_string.rstrip().split('\t')
    return bed[0], int(bed[1]), int(bed[2]), bed[3], bed[4], bed[5]


def annotate_bed_multi_core(cores, lines, stranded, transcript_priority,
                            gene_priority, cds_dict, features_dict,
                            transcript_id_key, type_key, out_file):
    """
    Instead of a single line, annotate on all lines using a
    multiprocessing module.

    :param cores: number of cores to use
    :param lines: list[string]
        list of lines from a BED file
    :param stranded: boolean
        if True, search for features on the same strand only. If
        False, search for same-stranded features, but also search
        opposite-stranded features if there is no overlapping
        same-strand feature available. If strand isn't specified
        in the BED file, search positive stranded features first.
    :param transcript_priority: list
        list of lists [['protein_coding','cds'],...] which explicitly
        order the priority by which features are annotated among
        transcripts.
    :param gene_priority: list
        list of lists [['protein_coding','cds'],...] which explicitly
        order the priority by which features are annotated among
        genes.
    :param cds_dict: dict
        dictionary of CDS start/end positions for each transcript
        see: get_all_cds_dict()
    :param features_dict: dict
        multilevel dictionary containing all features in all specified
        chromosomes.
    :param transcript_id_key: string
        identifier in GTF file that specifies transcript_id
        (ie. transcript_id)
    :param type_key: string
        identifier in GTF file that specifies the biotype
        (ie. transcript_type)
    :param out_file: string
        output file to write to disk.
    :return:
    """
    ### create shared-in-memory objects
    manager = multiprocessing.Manager()
    shared_cds_dict = manager.dict(cds_dict)
    shared_features_dict = manager.dict(features_dict)

    ### create pool
    print("setting up pool to run on {} cores".format(cores))
    po = multiprocessing.Pool(cores)
    fs = []
    # with open(out_file, 'w') as o:
    o = open(out_file, 'w')
    # for interval in bed_tool:  # for each line in bed file

    fs = []
    o = open(out_file, 'w')
    progress = trange(len(lines), desc='adding to pool', leave=False)
    with concurrent.futures.ProcessPoolExecutor(cores) as executor:
        for line in lines:
            chrom, start, end, name, score, strand = split_bed_interval(line)
            fs.append(executor.submit(
                annotate, chrom, start, end, name, score, strand,
                stranded, transcript_priority, gene_priority,
                shared_features_dict, shared_cds_dict, transcript_id_key, type_key
            ))
            progress.update(1)
    progress = trange(len(lines), desc='writing to disk', leave=False)
    for f in concurrent.futures.as_completed(fs):
        o.write('\t'.join([str(r) for r in f.result()]) + "\n")
        progress.update(1)
    o.close()


def annotate_bed_single_core(line, stranded, transcript_priority,
                            gene_priority, cds_dict, features_dict,
                            transcript_id_key, type_key, progress):
    """
    Annotates a bed6-formatted line given dictionaries containing
    cds and genic positions, and priority lists.

    :param line: string

    :param stranded: boolean
        if True, search for features on the same strand only. If
        False, search for same-stranded features, but also search
        opposite-stranded features if there is no overlapping
        same-strand feature available. If strand isn't specified
        in the BED file, search positive stranded features first.
    :param transcript_priority: list
        list of lists [['protein_coding','cds'],...] which explicitly
        order the priority by which features are annotated among
        transcripts.
    :param gene_priority: list
        list of lists [['protein_coding','cds'],...] which explicitly
        order the priority by which features are annotated among
        genes.
    :param cds_dict: dict
        dictionary of CDS start/end positions for each transcript
        see: get_all_cds_dict()
    :param features_dict: dict
        multilevel dictionary containing all features in all specified
        chromosomes.
    :param transcript_id_key: string
        identifier in GTF file that specifies transcript_id
        (ie. transcript_id)
    :param type_key: string
        identifier in GTF file that specifies the biotype
        (ie. transcript_type)

    :param progress: tqdm.tqdm()
        tqdm progress object
    :return:
    """
    chrom, start, end, name, score, strand = split_bed_interval(line)
    result = annotate(
        chrom, start, end, name, score, strand,
        stranded, transcript_priority, gene_priority,
        features_dict, cds_dict,
        transcript_id_key, type_key
    )
    progress.update(1)
    return result

def annotate_bed(db_file, bed_file, out_file, stranded, chroms,
             transcript_priority, gene_priority, species, append_chr, fuzzy,
             cores):
    """
    Given a bed6 file, return the file with an extra column containing
    '|' delimited gene annotations

    :param db_file: basestring
        gffutils database sqlite file
    :param bed_file: basestring
        unannotated bed file
    :param out_file: basestring
        output file name
    :param chroms: list
        list of strings pertaining to chromosomes we want to use to annotate.
    :param species: basestring
        'hg19' or 'mm10' or 'ce11' to specify the GTF annotation syntax.
    :param append_chr: boolean
        True if we need to add 'chr' to the database annotations
        For example, wormbase/Ensembl annotations use I/II/1/2/etc.
        instead of chrI/chr1, so we need to let the database know that.
    :return:
    """
    geneid_to_name_dict, exons_dict, transcripts_dict, \
    cds_dict, features_dict, cds_key, utr3_key, utr5_key, utr_key,\
    gene_name_key, transcript_id_key, type_key = create_definitions(
        db_file, chroms=chroms, species=species, append_chr=append_chr, fuzzy=fuzzy
    )

    i = open(bed_file, 'r')
    lines = i.readlines()
    progress = trange(len(lines))
    with open(out_file, 'w') as o:
        writer = csv.writer(o, delimiter='\t')
        writer.writerows([
            annotate_bed_single_core(  # TODO: implement multicore method
                line, stranded, transcript_priority,
                gene_priority, cds_dict, features_dict,
                transcript_id_key, type_key, progress
            ) for line in lines]
        )
    i.close()
    return 0