#!/usr/bin/env python

# Calculates the average sizes of each region


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
import os
import pybedtools
import pandas as pd
from . import annotate_bed
from tqdm import trange
from collections import defaultdict


def calculate_avg_cds_length(db, keys):
    cds_feature_tx_lengths = defaultdict(list)
    for cds_feature in db.features_of_type(keys['cds']):
        for transcript in cds_feature.attributes['transcript_id']:
            cds_feature_tx_lengths[transcript].append(cds_feature.end - cds_feature.start)

    cds_features = []
    for transcript in cds_feature_tx_lengths.keys():
        cds_features.append(sum(cds_feature_tx_lengths[transcript]))
    print('cds lengths: {}'.format(sum(cds_features) / float(len(cds_features))))

def calculate_avg_utr_lengths(db, cds_dict, keys):
    # five_prime_utr_features = defaultdict(list)
    # three_prime_utr_features = defaultdict(list)
    five_prime_utr_tx_lengths = defaultdict(list)
    three_prime_utr_tx_lengths = defaultdict(list)

    # Use CDS dict to determine for each UTR, whether it's 3/5'
    for utr_feature in db.features_of_type(keys['utr']):
        classified_utr = annotate_bed.classify_utr(utr_feature, cds_dict)
        if classified_utr == '5utr':
            for transcript in utr_feature.attributes['transcript_id']:
                five_prime_utr_tx_lengths[transcript].append(utr_feature.end - utr_feature.start)
        elif classified_utr == '3utr':
            for transcript in utr_feature.attributes['transcript_id']:
                three_prime_utr_tx_lengths[transcript].append(utr_feature.end - utr_feature.start)
    # Calculate the total lengths for all transcripts
    # five_prime_utr_tx_lengths = defaultdict(list)
    # three_prime_utr_tx_lengths = defaultdict(list)
    # for transcript in five_prime_utr_features.keys():
    #     for i in five_prime_utr_features[transcript]:
    #         five_prime_utr_tx_lengths[transcript].append(i.end - i.start)
    # for transcript in three_prime_utr_features.keys():
    #     for i in three_prime_utr_features[transcript]:
    #         three_prime_utr_tx_lengths[transcript].append(i.end - i.start)

    five_prime_utr_lengths = []
    three_prime_utr_lengths = []

    for transcript in three_prime_utr_tx_lengths.keys():
        three_prime_utr_lengths.append(sum(three_prime_utr_tx_lengths[transcript]))
    for transcript in five_prime_utr_tx_lengths.keys():
        five_prime_utr_lengths.append(sum(five_prime_utr_tx_lengths[transcript]))

    print('three_prime_utr_lengths: {}'.format(
        sum(three_prime_utr_lengths) / float(len(three_prime_utr_lengths)))
    )
    print('five_prime_utr_lengths: {}'.format(
        sum(five_prime_utr_lengths) / float(len(five_prime_utr_lengths)))
    )

def calculate_avg_lengths(db_file, species):
    ### get the parsed regions ###

    keys = annotate_bed.get_keys(species)
    db = gffutils.FeatureDB(db_file)

    cds_dict = annotate_bed.get_all_cds_dict(db, keys['cds'])
    calculate_avg_cds_length(db, keys)

    calculate_avg_utr_lengths(db, cds_dict, keys)

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

if __name__ == "__main__":
    main()