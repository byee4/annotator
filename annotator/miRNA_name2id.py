#!/usr/bin/env python

# Converts a file that is labeled with miRNA names into one that is labeled
# with miRNA accessions

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
import pandas as pd
from collections import defaultdict


def get_precursor_id_from_namerow(row, name2id, name_col):
    """ Uses name2id dict to return the id with the associated name. """
    try:
        return '|'.join(name2id[row[name_col]])
    except KeyError:
        print('Cant find {} in the database!'.format(row[name_col]))
        return ''


def get_mature_id_from_idrow(row, precursor2mature, name_col):
    """ Uses name2id dict to return the id with the associated name. """
    try:
        return '|'.join(precursor2mature[row[name_col]])
    except KeyError:
        print('Cant find {} in the database!'.format(row[name_col]))
        return ''


def add_precursor_id_to_dataframe(df, name2id_dict, name_col):
    """ adds id values to pandas dataframe """
    assert 'miRNA id' not in df.columns
    df['miRNA id'] = df.apply(get_precursor_id_from_namerow, args=[name2id_dict, name_col], axis=1)
    return df


def add_mature_id_to_dataframe(df, precursor2mature_dict, name_col):
    print(df.head())
    df['mature'] = df.apply(get_mature_id_from_idrow, args=[precursor2mature_dict, name_col], axis=1)
    return df


def build_precursor_id_to_mature_id_dict(gffdb_file, custom_file):
    """
    Uses a gffdb file and a custom file to
    :param gffdb_file:
    :param custom_file:
    :return:
    """
    precursor2mature = defaultdict(list)
    if gffdb_file is not None:
        try:
            db = gffutils.FeatureDB(gffdb_file)
            for transcript in db.features_of_type('miRNA'):
                for mature_id in transcript.attributes['ID']:
                    for parent in transcript.attributes['Derives_from']:
                        precursor2mature[parent].append(mature_id)
        except Exception as e:
            print(e)
            print("Could not use the gffdb file to annotate mature ids!")
    if custom_file is not None:
        try:
            with open(custom_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):  # comment chars
                        id, _, mature_ids = line.rstrip('\n').split('\t')
                        mature_ids = mature_ids.split('|')
                        if id in precursor2mature.keys():
                            precursor2mature[id].append(mature_ids)
                        else:
                            precursor2mature[id] = mature_ids
        except Exception as e:
            print(e)
            print("Could not use the custom file to annotate mature ids!")

    return precursor2mature


def build_name_to_precursor_id_dict(gffdb_file, custom_file):
    """
    Returns {name:id} dictionary. 
    Names with more than 1 associated ID will be appended, delimited by |.
    
    :param gffdb_file: string
        path of the gffdb.
        Database must contain:
        'miRNA_primary_transcript' in featuretype field
        'Name' in feature.attributes.keys()
        'ID' in feature.attributes.keys()
        
    :return name2id_dict: dict
        {name:id}
    """
    # delim = '|'
    precursor_name2id = defaultdict(list)
    if gffdb_file is not None:
        try:
            db = gffutils.FeatureDB(gffdb_file)
            for transcript in db.features_of_type('miRNA_primary_transcript'):
                for name in transcript.attributes['Name']:
                    precursor_name2id[name].append(transcript.id)

        except Exception as e:
            print(e)
            print("Could not use the gffdb file to annotate precursor ids!")
    if custom_file is not None:
        try:
            with open(custom_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):  # comment chars
                        id, name = line.rstrip('\n').split('\t')[:2]  # ignore other columns
                        precursor_name2id[name].append(id)
        except Exception as e:
            print(e)
            print("Could not use the custom file to annotate precursor ids!")
    return precursor_name2id


def read_and_append_id(in_file, sep, name_col, gffdb_file, custom_file, add_mature, delete_original_column):
    """
    Reads in a delimited file with miRNA names in at least one column
    Returns dataframe containing matching miRNA precursor ids.
    
    :param in_file: string
        input delimited file
    :param sep: string
        delimiter
    :param name_col: string
        column name where miRNA names exist
    :param gffdb_file: string
        path to gffutils database file
    :param delete_original_column: bool
        if True, delete the original name column (specified with name_col)
        if False, keep it.
    :return: 
    """
    df = pd.read_table(in_file, sep=sep)

    precursor_name2id_dict = build_name_to_precursor_id_dict(gffdb_file, custom_file)
    df = add_precursor_id_to_dataframe(df, precursor_name2id_dict, name_col)

    if add_mature:
        precursor2mature_dict = build_precursor_id_to_mature_id_dict(gffdb_file, custom_file)
        df = add_mature_id_to_dataframe(df, precursor2mature_dict, 'miRNA id')

    if delete_original_column:
        del df[name_col]
        df.set_index('miRNA id', inplace=True)
    else:
        df.set_index(name_col, inplace=True)
    return df


def convert(in_file, sep, name_col, out_file, gtfdb_file, custom_file, add_mature, delete_original_column=True):
    df = read_and_append_id(
        in_file, sep, name_col, gtfdb_file, custom_file,
        add_mature,
        delete_original_column
    )

    df.to_csv(out_file, sep=sep, index=True, header=True)


def main():
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "--output",
        dest="output",
        help="output file",
        required=False
    )
    parser.add_argument(
        "--input",
        dest="input",
        help="input file",
        required=True
    )
    parser.add_argument(
        "--sep",
        dest="sep",
        help="default is tabbed (\'\\t\'), but can be specified here",
        required=False,
        default='\t'
    )
    parser.add_argument(
        "--name_col",
        dest="name_col",
        help="name of the column that contains the miRNA name",
        required=False,
        default=None
    )
    parser.add_argument(
        "--add_mature",
        action='store_true',
        help="add the mature sequences, | delimited, as another column",
        required=False,
        default=False
    )
    parser.add_argument(
        "--gffdb",
        dest="gffdb",
        help="gff database file create from gffutils (Currently only "
             "GFF files are allowed)",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--custom",
        dest="custom",
        help="custom tabbed file (id\tname) for mirs not in the gff database",
        required=False,
        default=None
    )
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    in_file = args.input
    sep = args.sep
    name_col = args.name_col
    out_file = args.output
    gffdb_file = args.gffdb
    custom_file = args.custom
    add_mature = args.add_mature

    convert(in_file, sep, name_col, out_file, gffdb_file, custom_file, add_mature)
if __name__ == "__main__":
    main()