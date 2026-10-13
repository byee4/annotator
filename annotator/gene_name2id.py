#!/usr/bin/env python

# Converts a file that is labeled with gene names into one that is labeled
# with gene ids

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


def get_id_from_namerow(row, name2id, name_col):
    """Look up ordered gene IDs for one table row. The lookup value is converted to str.

    :param row: Mapping or pandas row holding the lookup column.
    :param name2id: Mapping of names to ordered ID lists.
    :param name_col: Column name used for lookup.
    :return: A pipe-joined ID string, or an empty string for an unknown key.
    """
    try:
        return '|'.join(name2id[str(row[name_col])])
    except KeyError:
        print('Cant find {} in the database!'.format(str(row[name_col])))
        return ''

def add_id_to_dataframe(df, name2id_dict, name_col):
    """Append the 'gene id' lookup column to a DataFrame. Raises AssertionError if the ID column already exists.

    :param df: Pandas DataFrame to update in place.
    :param name2id_dict: Mapping of names to ordered ID lists.
    :param name_col: Column name used for lookup.
    :return: The same DataFrame with joined IDs; unknown keys have empty strings.
    """
    assert 'gene id' not in df.columns
    df['gene id'] = df.apply(get_id_from_namerow, args=[name2id_dict, name_col], axis=1)
    return df


def build_name_to_id_dict(gffdb_file, custom_file):
    """Build an identifier mapping from optional database and custom table sources.

    :param gffdb_file: Optional path to a gffutils database, or None.
    :param custom_file: Optional tab-delimited mapping path; comment lines start with #. Gene mappings use ID/name columns; microRNA mappings use precursor ID/name/pipe-delimited mature IDs.
    :return: A defaultdict(list) mapping names to ID lists; empty when both sources are None. Source errors are printed and do not propagate.
    """
    # delim = '|'
    name2id = defaultdict(list)
    if gffdb_file is not None:
        try:
            db = gffutils.FeatureDB(gffdb_file)
            for transcript in db.features_of_type('gene'):
                for gene_name in transcript.attributes['gene_name']:
                    for gene_id in transcript.attributes['gene_id']:
                        name2id[gene_name].append(gene_id)

        except Exception as e:
            print(e)
            print("Could not use the gffdb file to annotate precursor ids!")
    if custom_file is not None:
        try:
            with open(custom_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):  # comment chars
                        gene_id, gene_name = line.rstrip('\n').split('\t')  # ignore other columns
                        name2id[gene_name].append(gene_id)
        except Exception as e:
            print(e)
            print("Could not use the custom file to annotate precursor ids!")
    return name2id


def read_and_append_id(in_file, sep, name_col, gffdb_file, custom_file, delete_original_column):
    """Read a header-bearing table and append gene IDs.

    :param in_file: Input table path with a header row.
    :param sep: Input/output table delimiter.
    :param name_col: Column name used for lookup.
    :param gffdb_file: Optional path to a gffutils database, or None.
    :param custom_file: Optional tab-delimited mapping path; comment lines start with #. Gene mappings use ID/name columns; microRNA mappings use precursor ID/name/pipe-delimited mature IDs.
    :param delete_original_column: True to drop the name column and index by ID; False to index by the original name.
    :return: A DataFrame indexed by original names or IDs according to delete_original_column.
    """
    df = pd.read_table(in_file, sep=sep)

    name2id_dict = build_name_to_id_dict(gffdb_file, custom_file)
    df = add_id_to_dataframe(df, name2id_dict, name_col)

    if delete_original_column:
        del df[name_col]
        df.set_index('gene id', inplace=True)
    else:
        df.set_index(name_col, inplace=True)
    return df


def convert(in_file, sep, name_col, out_file, gtfdb_file, custom_file, delete_original_column):
    """Convert a header-bearing identifier table and write its indexed result.

    :param in_file: Input table path with a header row.
    :param sep: Input/output table delimiter.
    :param name_col: Column name used for lookup.
    :param out_file: Output table path; includes header and index.
    :param gtfdb_file: Optional path to a gffutils database, or None.
    :param custom_file: Optional tab-delimited mapping path; comment lines start with #. Gene mappings use ID/name columns; microRNA mappings use precursor ID/name/pipe-delimited mature IDs.
    :param delete_original_column: True to drop the name column and index by ID; False to index by the original name.
    :return: None; the output file contains the converted table.
    """
    df = read_and_append_id(
        in_file, sep, name_col, gtfdb_file, custom_file,
        delete_original_column
    )
    df.to_csv(out_file, sep=str(sep), index=True, header=True)


def main():
    # Setup argument parser
    """Convert gene names using command-line input/output paths and mapping options.

    Reads sys.argv: custom and/or gffdb supply ID mappings; name_col identifies
    the name column; sep defaults to tab. delete_original_column indexes by ID
    instead of name. Returns None after writing the converted table.
    """
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
        help="delimiter, default is tabbed (\'\\t\'), but can be specified here",
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
        help="custom tabbed file (id\tname) for genes not in the gff database",
        required=False,
        default=None
    )
    parser.add_argument(
        "--delete_original_column",
        dest="delete_original_column",
        help="if True, delete the original column (name_col)",
        required=False,
        action='store_true',
        default=False
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
    delete_original_column = args.delete_original_column

    convert(in_file, sep, name_col, out_file, gffdb_file, custom_file, delete_original_column)
if __name__ == "__main__":
    main()