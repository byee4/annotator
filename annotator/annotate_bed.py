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
import pybedtools
import cProfile
from . import annotation_functions as af

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


def create_definitions(db_file, chroms=[], gtf_format='gencode',
                       append_chr=False, fuzzy=0):
    """
    This method creates definitions and populates dictionaries to be used
    for annotation. Previously initialized in Annotator().

    :param db_file:
    :param chroms:
    :param gtf_format: basestring
    'gencode', 'ensembl', or 'refseq' to specify the GTF annotation syntax.
    :param append_chr:
    :param fuzzy:
    :return:
    """

    keys = af.get_keys(gtf_format)

    db = gffutils.FeatureDB(db_file)
    featuretypes = [x.lower() for x in list(db.featuretypes())]

    progress = trange(4, desc='Initializing/creating defs')
    progress.set_description("Building gene ids -> gene names dictionary")

    progress.update(1)
    if len(chroms) != 0:  # if use specific chromosomes, otherwise hash all
        chromosomes = set(chroms)
    else:
        chromosomes = af.chromosome_set(db)

    exons_dict = {}  # holds all of the transcript_id:exon start/stop
    transcripts_dict = {}  # holds all of the transcript_id:start/stop
    cds_dict = {}  # holds all of the transcript_id:cds start/stop

    # If 'exon' features exist in db, catalog it.
    if 'exon' in featuretypes:
        progress.set_description("Adding all exon boundaries")
        exons_dict = af.get_all_exons_dict(
            db, keys['exon'], keys['transcript_id']
        )
    progress.update(1)
    # If 'transcript' features exist in db, catalog it.
    if 'transcript' in featuretypes:
        progress.set_description("Adding all transcript boundaries")
        transcripts_dict = af.get_all_transcripts_dict(
            db, keys['transcript'], keys['transcript_id']
        )
    progress.update(1)
    # If 'cds' features exist in db, catalog it.
    if 'cds' in featuretypes:
        progress.set_description("Adding all CDS boundaries")
        cds_dict = af.get_all_cds_dict(db, keys['cds'])
    progress.update(1)
    progress.set_description("Adding all features in DB to dictionary")
    features_dict, num_features = hash_features(
        db, chromosomes, append_chr, fuzzy,
        exons_dict, transcripts_dict,
        keys['transcript_id'], featuretypes
    )

    return exons_dict, transcripts_dict, cds_dict, features_dict, keys


def hash_features(db, chromosomes, append_chr, fuzzy,
                  exons_dict, transcripts_dict,
                  transcript_id_key, featuretypes):
    """
    hashes features by position.
    :return chr19_features_dict : collections.defaultdict()
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

            # TODO: implement 'fuzzy gene matching' feature
            # if element.featuretype == 'gene':
            #     if element.strand == '+':
            #         element.start = element.start - fuzzy
            #     elif element.strand == '-':
            #         element.end = element.end + fuzzy
            # if element.featuretype == 'gene':
            #     if element.strand == '+':
            #         element.end = element.end + fuzzy
            #     elif element.strand == '-':
            #         element.start = element.start - fuzzy

            introns = []
            # if we see at least one intron, we assume introns are annotated!
            if element.featuretype == 'transcript' and \
                            'intron' not in featuretypes:
                for transcript_id in element.attributes[transcript_id_key]:
                    exons = exons_dict[transcript_id]
                    transcript = transcripts_dict[transcript_id]
                    # Find introns
                    introns = af.find_introns(transcript, exons)

            # Create multi-key hash and store elements
            start = int(element.start / HASH_VAL)
            end = int(element.end / HASH_VAL)

            for i in range(start, end + 1):
                element.chrom = fixed_chrom
                append = features_dict[fixed_chrom, i, element.strand].append

                # append elements to features dictionary
                append(element)

                # append inferred introns to features dictionary
                for intron in introns:
                    #  TODO: figure out a way to remove deepcopy requirement
                    intron_interval = pybedtools.create_interval_from_list(
                        [element.chrom, str(int(intron['start'])-1), intron['end'],  # turn 1 based into 0 based
                        'intron', str(element.score), element.strand]
                    )
                    proxdist_dict = af.get_proxdist_from_intron(intron_interval)
                    for prox in proxdist_dict['prox']:
                        proxintron_feature = copy.deepcopy(element)
                        proxintron_feature.start = prox.start + 1
                        proxintron_feature.end = prox.end
                        proxintron_feature.featuretype = prox.name
                        append(proxintron_feature)
                    for dist in proxdist_dict['dist']:
                        distintron_feature = copy.deepcopy(element)
                        distintron_feature.start = dist.start + 1
                        distintron_feature.end = dist.end
                        distintron_feature.featuretype = dist.name
                        append(distintron_feature)

            num_features += 1
        progress.update(1)
    # for feature in chr19_features_dict['chr7', 129, '+']:
    #     if feature.attributes['gene_id'][0] == 'ENSG00000128607.9':
    #         if feature.featuretype == 'CDS':
    #             print(feature)
    return features_dict, num_features


def get_all_overlapping_features_from_query_stranded(
        chrom, qstart, qend, strand, features_dict):
    """
    Given a query location (chr, start, end), return all features that
    overlap by at least one base. Functions similarly to gffutils db.region(),
    but uses the pre-hashed self.chr19_features_dict to greatly speed things up.

    :param chrom : string
    :param qstart : int
    :param qend : int
    :param strand : string
    :param features_dict : dict
        dictionary of positions for all features
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
            # query completely contains feature
            if qstart <= feature.start and qend >= feature.end:
                feature.attributes['overlap'] = 'query_contains_feature'
                append(feature)
            # feature completely contains query
            elif qstart >= feature.start and qend <= feature.end:
                feature.attributes['overlap'] = 'feature_contains_query'
                append(feature)
            # feature partially overlaps (qstart < fstart < qend)
            elif qstart <= feature.start and qend >= feature.start:
                feature.attributes[
                    'overlap'] = 'partial_by_{}_bases'.format(
                    qend - (feature.start - 1))
                append(feature)
            # feature partially overlaps (qstart < fend < qend)
            elif qstart <= feature.end and qend >= feature.end:
                feature.attributes[
                    'overlap'] = 'partial_by_{}_bases'.format(
                    feature.end - qstart)
                append(feature)

    return features


def get_all_overlapping_features_from_query_unstranded(
        chrom, qstart, qend, strand, features_dict):
    """
    Returns overlapping features from a query.

    :param chrom: string
        query chromosome/source
    :param qstart: int
        query start
    :param qend: int
        query end
    :param strand: string
        query strand
    :param features_dict: dict
        dictionary of positions for all features
    :return:
    """
    if strand is None:  # If no strand given, return features from both strands
        positive_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, '+', features_dict
        )
        negative_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, '-', features_dict
        )
        return positive_features + negative_features
    # return features on the same strand, only return features on
    # opposite strand if it can't find anything otherwise.
    else:

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
    If stranded is True, then return only the features that overlap on
    the same strand. Otherwise, features returned need not to be with respect
    to strand (although if strand is given, it will return same-stranded
    features first before attempting to find those on the opposite strand).

    :param chrom:
    :param qstart:
    :param qend:
    :param strand:
    :param features_dict:
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
            # TODO: implement better. For now, returns both strands
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

    :param transcript_type: string
    :return is_protein_coding: boolean

    """
    if transcript_type == 'protein_coding':
        return True

    return False


def classify_transcript_type(
        transcript_type, gene_type, distinct_nc_tx=None
):
    """
    deprecates is_protein_coding().

    :param transcript_type: string
    :param gene_type: string
    :param distinct_nc_tx: list or None
        We can pull out any distinct noncoding regions using this parameter.
        For example, we can specify distinct_nc_tx = ['lincRNA', 'pseudogene'],
        and this function will NOT bin these categories into the catchall
        'non_coding' category.
    :return transcript_type_class: string

    """

    # if the transcript type is protein coding, easy. return 'protein coding'.
    if transcript_type == 'protein_coding':
        return 'protein_coding'
    # TODO: Is this protein coding?
    # elif gene_type == 'protein_coding' and transcript_type == 'retained_intron':
    #     return 'protein_coding'
    # if we specify any non_coding type to pull out, returns the type here.
    # if transcript_type == 'snoRNA':
    #     return 'snoRNA'
    # if transcript_type == 'scaRNA':
    #     return 'scaRNA'
    if transcript_type == 'miRNA':
        return 'miRNA'
    if distinct_nc_tx is not None:
        for distinct_type in distinct_nc_tx:
            if transcript_type == distinct_type:
                return distinct_type
    # if not either protein_coding or anything in distinct_nc_tx, return non_coding.
    return 'non_coding'


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
        # print("FEATURE STRING: {}".format(feature_string))
        transcript, start, end, strand, feature_type, gene_id, \
        gene_name, transcript_type_list, gene_type_list, overlap = feature_string.split(':')
        transcript_type_list = transcript_type_list.split(',')
        for transcript_type, gene_type in zip(transcript_type_list, gene_type_list):
            combined_dict[
                classify_transcript_type(transcript_type, gene_type),
                feature_type
            ].append(
                feature_string
            )
    # print("COMBINED DICT STEP 1: {}".format(combined_dict))
    # return the highest one
    combined_dict = OrderedDict(
        combined_dict
    )  # turn into ordered dict, is that ok?
    # print("combined dict: {}".format(combined_dict))
    # print("COMBINED DICT STEP 2: {}".format(combined_dict))
    try:
        # append nonexisting regions to the end of the priority list
        for key, _ in iteritems(combined_dict):
            if list(key) not in priority:
                priority.append(list(key))
        # print("PRIORITY: ",priority)
        combined_dict = sorted(  # sort based on priority list
            iteritems(combined_dict),
            key=lambda x: priority.index([x[0][0], x[0][1]])
        )
        # print("COMBINED DICT STEP 3: {}".format(combined_dict))
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
        dictionary where chr19_keys = transcript ids and values = list of
        features associated with each transcript that overlaps with a query
        (ie. {'ENSTXYZ':['CDS', 'exon', 'transcript']}
    :param transcript_priority: list[tuple]
        list of tuples (featuretype, transcripttype)
    :return unique_genes: defaultdict[list]
    """
    unique_genes = defaultdict(list)
    # print("unique transcript features: {}".format(unique_transcript_features))
    if len(unique_transcript_features.keys()) == 0:
        return 'intergenic'
    # print('all transcript features: {}'.format(unique_transcript_features.keys()))
    for transcript in unique_transcript_features.keys():
        top_transcript = return_highest_priority_feature(
            unique_transcript_features[transcript],
            transcript_priority
        )[1][0]  # [0] contains the dictionary key
        gene_list = top_transcript.split(':')[5].split(',')
        for gene in gene_list:
            unique_genes[gene].append(top_transcript)

    for g in unique_genes.keys():
        unique_genes[g] = return_highest_priority_feature(
            unique_genes[g],
            transcript_priority
        )[1]
    return unique_genes


def prioritize_genes(unique_genes, gene_priority):
    """
    Given a dictionary of genes and gene type priority, return 
    the highest priority gene whose feature type is the highest priority.
    
    :param unique_genes: 
    :param gene_priority: 
    :return: 
    """
    # print("UNIQUE GENES", unique_genes)
    final_transcripts = []
    final_transcripts_append = final_transcripts.append
    for gene, transcripts in iteritems(unique_genes):
        for transcript in transcripts:
            final_transcripts_append(transcript)

    feature_type, final_genes = return_highest_priority_feature(
        final_transcripts, gene_priority
    )
    if feature_type[0] == 'miRNA':  # TODO: fix. hacky
        ff = []
        for f in sorted(final_genes):
            f = f.replace(
                'exon', 'miRNA'
            )
            ff.append(f)
        return ff
    elif feature_type[0] == 'non_coding':  # TODO: fix. hacky
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

    # Filter gene features (we deal with priority on a per-transcript basis)
    for feature_string in formatted_features:
        if feature_string.split(':')[4] != parent:
            transcript = feature_string.split(':')[0]
            unique_transcripts[transcript].append(
                feature_string
            )
    # Get one transcript per gene whose featuretype has highest priority.
    unique_genes = prioritize_transcript(
        unique_transcripts, transcript_priority
    )
    # Get one gene whose featuretype has highest priority
    top_gene = prioritize_genes(unique_genes, gene_priority)
    return top_gene


def split_string(priority_string, delimiter=':'):
    """
    Parses/splits a string with an expected format and returns
    three important fields: region, gene_id, gene_name. Returns
    'intergenic' for all three fields if intergenic

    :param priority_string: string
        colon-delimited string containing region geneid and
        genename in fields 4-6 (0-based open).
    :param delimiter: change
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
        features_dict, cds_dict, keys):

    """
    Given position parameters, return the annotation string.

    :param chrom:
    :param start:
    :param end:
    :param name:
    :param score:
    :param strand:
    :param stranded:
    :param transcript_priority:
    :param gene_priority:
    :param features_dict:
    :param cds_dict:
    :param keys:
    :return:
    """
    append_count = 0  # number of features that were found to overlap reigon

    overlapping_features = get_all_overlapping_features_from_query(
        chrom, start, end, strand, features_dict, stranded
    )
    to_append = ''  # full list of genes overlapping features
    transcript = defaultdict(list)
    for feature in overlapping_features:  # for each overlapping feature
        # for M10 annotations, I see gene featuretypes without
        # transcript_ids. That's fine. We can re-construct genes from the
        # transcript features, so ignore anything that doesn't have a
        # transcript_id in it.
        if keys['transcript_id'] in feature.attributes.keys():
            for transcript_id in feature.attributes[
                keys['transcript_id']
            ]:  # multiple genes can be associated with one feature
                transcript[transcript_id].append(
                    feature)  # append features to their respective genes

    # If we find no overlapping features, return 'intergenic'
    if len(transcript.keys()) == 0:
        return chrom, start, end, name, score, strand, \
               'intergenic', 'intergenic', \
               'intergenic', 'intergenic'
    # Generate an annotation string
    for transcript, features in iteritems(transcript):
        for feature in features:
            append_count += 1
            # if we only have UTR, specify what kind of UTR (3' or 5').
            if feature.featuretype is None:
                print("Warning: featuretype for feature {} is None".format(feature))
            elif feature.featuretype == keys['utr']:
                feature.featuretype = af.classify_utr(feature, cds_dict)
            to_append += "{}:{}:{}:{}:{}:".format(
                transcript,
                feature.start - 1,  # report 0 based
                feature.end,
                feature.strand,
                feature.featuretype,
            )
            if keys['gene_id'] in feature.attributes.keys():
                for t in feature.attributes[keys['gene_id']]:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            else:
                to_append = to_append + '-:'
            if keys['gene_name'] in feature.attributes.keys():
                for t in feature.attributes[keys['gene_name']]:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            elif keys['transcript_id'] in feature.attributes.keys():
                for t in feature.attributes[keys['transcript_id']]:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            else:
                to_append = to_append + '-:'
            if keys['transcript_type'] in feature.attributes.keys():
                for t in feature.attributes[keys['transcript_type']]:
                    to_append += '{},'.format(t)
                to_append = to_append[:-1] + ':'
            else:
                to_append = to_append + '-:'
            if keys['gene_type'] in feature.attributes.keys():
                for t in feature.attributes[keys['gene_type']]:
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

    to_append = to_append[:-1]  # remove the trailing ':'
    if append_count == 1:  # if just one region, just return it.
        region, gene, rname = split_string(to_append)
    else:
        priority = prioritize_transcript_then_gene(
            parse_annotation_string(to_append),
            transcript_priority,
            gene_priority
        )
        region, gene, rname = split_string(priority[0])
        for highest_priority_gene in priority[1:]:
            current_region, current_gene, current_rname = split_string(
                highest_priority_gene
            )
            # there should really just be one highest priority region.
            assert region == current_region
            # if the priority region contains multiple genes, return all
            if gene != current_gene:
                gene = gene + ',' + current_gene
            if rname != current_rname:
                rname = rname + ',' + current_rname
    return chrom, start, end, name, score, strand, \
           gene , rname, region, to_append


def split_bed_interval(bed_string):
    """
    Takes a bed file as a line (tab delimited) and returns 6 fields
    correctly typed.
    :param bed_string: string
    :return chrom, start, end, name, score, strand: string
        bed interval
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
    # create shared-in-memory objects
    manager = multiprocessing.Manager()
    shared_cds_dict = manager.dict(cds_dict)
    shared_features_dict = manager.dict(features_dict)

    # create pool
    print("setting up pool to run on {} cores".format(cores))

    fs = []
    o = open(out_file, 'w')
    progress = trange(len(lines), desc='adding to pool', leave=False)
    with concurrent.futures.ProcessPoolExecutor(cores) as executor:
        for line in lines:
            chrom, start, end, name, score, strand = split_bed_interval(line)
            fs.append(executor.submit(
                annotate, chrom, start, end, name, score, strand,
                stranded, transcript_priority, gene_priority,
                shared_features_dict, shared_cds_dict, transcript_id_key,
                type_key
            ))
            progress.update(1)
    progress = trange(len(lines), desc='writing to disk', leave=False)
    for f in concurrent.futures.as_completed(fs):
        o.write('\t'.join([str(r) for r in f.result()]) + "\n")
        progress.update(1)
    o.close()


def annotate_bed_single_core(line, stranded, transcript_priority,
                             gene_priority, cds_dict, features_dict,
                             keys, progress):
    """
    Annotates a bed6-formatted line given dictionaries containing
    cds and genic positions, and priority lists.

    :param line: string
        a single BED6-style line
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
    :param keys: dict
        dictionary of chr19_keys that specify the GTF language.
    :param progress: tqdm.tqdm()
        tqdm progress object
    :return result: string
        The resultant BED file, plus all the annotation columns
    """
    chrom, start, end, name, score, strand = split_bed_interval(line)
    result = annotate(
        chrom, start, end, name, score, strand,
        stranded, transcript_priority, gene_priority,
        features_dict, cds_dict, keys
    )
    progress.update(1)
    return result


def annotate_bed(db_file, bed_files, out_files, stranded, chroms,
             transcript_priority, gene_priority, gtf_format, append_chr, fuzzy,
             cores):
    """
    Given a bed6 file, return the file with an extra column containing
    '|' delimited gene annotations

    :param db_file: basestring
        gffutils database sqlite file
    :param bed_file: list
        list of unannotated bed files
    :param out_file: list
        list of output file names
    :param chroms: list
        list of strings pertaining to chromosomes we want to use to annotate.
    :param gtf_format: basestring
        'gencode', 'ensembl', or 'refseq' to specify the GTF annotation syntax.
    :param append_chr: boolean
        True if we need to add 'chr' to the database annotations
        For example, wormbase/Ensembl annotations use I/II/1/2/etc.
        instead of chrI/chr1, so we need to let the database know that.
    :return:
    """

    exons_dict, transcripts_dict, \
    cds_dict, features_dict, keys = create_definitions(
        db_file, chroms=chroms, gtf_format=gtf_format, append_chr=append_chr,
        fuzzy=fuzzy
    )
    for n in range(len(bed_files)):
        bed_file = bed_files[n]
        out_file = out_files[n]
        i = open(bed_file, 'r')
        lines = i.readlines()
        progress = trange(len(lines))
        with open(out_file, 'w') as o:
            writer = csv.writer(o, delimiter='\t')
            writer.writerows([
                annotate_bed_single_core(  # TODO: implement multicore method
                    line, stranded, transcript_priority,
                    gene_priority, cds_dict, features_dict,
                    keys, progress
                ) for line in lines]
            )
        i.close()
    return 0
