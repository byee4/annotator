#!/usr/bin/env python

# Calculates the average and total sizes of each region


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
from . import annotation_functions as af
from . import create_region_bedfiles
from collections import defaultdict


def calculate_total_cds_length(db, keys):
    """
    Calculates the total length of all CDS features in a db
    
    :param db: gffutils.FeatureDB
    :param keys: dict
        {'region':region defined in GTF}. Sets the naming scheme 
         for GTF parsing
    :return total_coverage: int
    """

    ### This method counts every CDS feature, may double count overlapping features ###

    # total_cds_feature_tx_lengths = 0
    # for cds_feature in db.features_of_type(keys['cds']):
    #     total_cds_feature_tx_lengths += (cds_feature.end - cds_feature.start)
    # print('cds lengths (total): {}'.format(total_cds_feature_tx_lengths))

    ### This method counts only nonoverlapping CDS regions ###
    cds = create_region_bedfiles.create_cds_region_bedfile(db, keys, None, True)
    cds = cds.sort()
    print('cds lengths (total): {}'.format(cds.total_coverage()))
    return cds.total_coverage()


def calculate_avg_cds_length(db, keys):
    """
    Calculates the average length of CDS regions for each transcript.
    
    :param db: gffutils.FeatureDB
    :param keys: dict
        {'region':region defined in GTF}. Sets the naming scheme 
         for GTF parsing
    :return avg_cds_len: float
    """

    cds_feature_tx_lengths = defaultdict(list)

    # get all cds regions that belong to each transcript
    for cds_feature in db.features_of_type(keys['cds']):
        for transcript in cds_feature.attributes['transcript_id']:
            cds_feature_tx_lengths[transcript].append(cds_feature.end - cds_feature.start)

    # calculate for each transcript, the total cds length
    cds_features = []
    for transcript in cds_feature_tx_lengths.keys():
        cds_features.append(sum(cds_feature_tx_lengths[transcript]))
    print('cds lengths (avg): {}'.format(sum(cds_features) / float(len(cds_features))))
    return sum(cds_features) / float(len(cds_features))


def calculate_total_utr_lengths(db, cds_dict, keys):
    """
    Calculates the total length of all UTR features in a db

    :param db: gffutils.FeatureDB
    :param cds_dict: collections.defaultdict
        defaultdict{transcript:{'start':START, 'end':END}}
    :param keys: dict
        {'region':region defined in GTF}. Sets the naming scheme 
         for GTF parsing
    :return total_coverage: int
    """

    total_five_prime_utr_tx_lengths = 0
    total_three_prime_utr_tx_lengths = 0

    ### This method counts every UTR feature - may double count overlapping features ###
    # Use CDS dict to determine for each UTR, whether it's 3/5'
    # for utr_feature in db.features_of_type(keys['utr']):
    #     classified_utr = annotate_bed.classify_utr(utr_feature, cds_dict)
    #     if classified_utr == '5utr':
    #         total_five_prime_utr_tx_lengths += (utr_feature.end - utr_feature.start)
    #     elif classified_utr == '3utr':
    #         total_three_prime_utr_tx_lengths += (utr_feature.end - utr_feature.start)


    ### This method counts only nonoverlapping UTR regions ###
    utr5_out = None  # no need to write to file
    utr3_out = None  # no need to write to file
    by_transcript = True  # dont merge per gene, keep transcript separate

    utr5, utr3 = create_region_bedfiles.create_utr_region_bedfiles(
        db, keys, cds_dict, utr3_out, utr5_out, by_transcript
    )
    utr5 = utr5.sort()
    utr3 = utr3.sort()

    print('five_prime_utr_lengths (total): {}'.format(utr5.total_coverage()))
    print('three_prime_utr_lengths (total): {}'.format(utr3.total_coverage()))
    return total_five_prime_utr_tx_lengths, total_three_prime_utr_tx_lengths


def calculate_avg_utr_lengths(db, cds_dict, keys):
    """
    Calculates the average length of all UTR features in a db
        
    :param db: gffutils.FeatureDB
    :param cds_dict: collections.defaultdict
        defaultdict{transcript:{'start':START, 'end':END}}
    :param keys: dict
        {'region':region defined in GTF}. Sets the naming scheme 
         for GTF parsing 
    :return: 
    """

    # {transcript_id: [len(tx1), len(tx2), ...]
    five_prime_utr_tx_lengths = defaultdict(list)
    three_prime_utr_tx_lengths = defaultdict(list)

    # Use CDS dict to determine length for all UTRs, whether it's 3/5'
    for utr_feature in db.features_of_type(keys['utr']):
        classified_utr = af.classify_utr(utr_feature, cds_dict)
        if classified_utr == '5utr':
            for transcript in utr_feature.attributes['transcript_id']:
                five_prime_utr_tx_lengths[transcript].append(
                    utr_feature.end - utr_feature.start
                )
        elif classified_utr == '3utr':
            for transcript in utr_feature.attributes['transcript_id']:
                three_prime_utr_tx_lengths[transcript].append(
                    utr_feature.end - utr_feature.start
                )

    # list of sums of each transcript utr length
    five_prime_utr_lengths = []
    three_prime_utr_lengths = []

    for transcript in three_prime_utr_tx_lengths.keys():
        three_prime_utr_lengths.append(
            sum(three_prime_utr_tx_lengths[transcript])
        )
    for transcript in five_prime_utr_tx_lengths.keys():
        five_prime_utr_lengths.append(
            sum(five_prime_utr_tx_lengths[transcript])
        )

    # divide by total length by number of transcripts
    avg_three_prime_utr_length = sum(three_prime_utr_lengths) / float(len(three_prime_utr_lengths))
    avg_five_prime_utr_length = sum(five_prime_utr_lengths) / float(len(five_prime_utr_lengths))

    print(
        'three_prime_utr_lengths (avg): {}'.format(
            avg_three_prime_utr_length
        )
    )
    print(
        'five_prime_utr_lengths (avg): {}'.format(
            avg_five_prime_utr_length
        )
    )


def calculate_total_intron_lengths(db, exons_dict, transcripts_dict, keys):
    """
    Calculates the total length of distal, proximal, and all introns.
    
    :param db: 
    :param exons_dict: 
    :param transcripts_dict: 
    :param keys: 
    :return: 
    """

    proxintron_out = None  # don't write to any file
    distintron_out = None  # don't write to any file
    allintron_out = None  # don't write to any file
    by_transcript = True  # don't merge by gene (keep transcripts separate)

    proxintrons, distintrons, allintrons = create_region_bedfiles.create_intron_region_bedfiles(
        db, exons_dict, transcripts_dict, keys,
        proxintron_out, distintron_out, allintron_out, by_transcript
    )
    proxintrons = proxintrons.sort()
    distintrons = distintrons.sort()
    allintrons = allintrons.sort()

    print('proximal intron lengths (total): {}'.format(proxintrons.total_coverage()))
    print('distal intron lengths (total): {}'.format(distintrons.total_coverage()))
    print('all intron lengths (total): {}'.format(allintrons.total_coverage()))


def calculate_avg_lengths(db_file, species):
    """
    Prints all average lengths for CDS and 3/5'UTRs.
    Can be used to generate ratios for metagene plots.
    
    :param db_file: gffutils.FeatureDB
    :param species: string
        used to generate keys corresonding to GTF column nomenclature.
    :return: 
    """

    keys = af.get_keys(species)
    db = gffutils.FeatureDB(db_file)

    cds_dict = af.get_all_cds_dict(db, keys['cds'])

    calculate_avg_cds_length(db, keys)
    calculate_avg_utr_lengths(db, cds_dict, keys)


def calculate_total_lengths(db_file, species):
    """"
    Prints number of bases labeled as CDS, UTR, and intron in a gtf db file.
    Can be used to generate ratios for region distribution pie charts.
    
    :param db_file: gffutils.FeatureDB
    :param species: string
        used to generate keys corresonding to GTF column nomenclature.
    :return: 
    """
    keys = af.get_keys(species)
    db = gffutils.FeatureDB(db_file)

    cds_dict = af.get_all_cds_dict(db, keys['cds'])
    exons_dict = af.get_all_exons_dict(db, keys['transcript_id'])
    transcripts_dict = af.get_all_transcripts_dict(db, keys[
        'transcript_id'])

    calculate_total_cds_length(db, keys)
    calculate_total_utr_lengths(db, cds_dict, keys)
    calculate_total_intron_lengths(db, exons_dict, transcripts_dict, keys)


def main():
    # Setup argument parser
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
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    db_file = args.db_file
    species = args.species

    calculate_avg_lengths(db_file, species)
    calculate_total_lengths(db_file, species)
if __name__ == "__main__":
    main()