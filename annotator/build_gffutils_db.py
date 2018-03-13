#!/usr/bin/env python

# Builds a FeatureDB database file using gffutils. For use with anything in
# annotator.

# these two are really a minimum

from __future__ import print_function
from __future__ import division

# uncomment from this compatibility import list, as py3/py2 support progresses

from __future__  import absolute_import
from __future__  import unicode_literals
from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
from future.utils import iteritems

from argparse import ArgumentParser
import sys

import gffutils

def build_db(
        annotation_file, db_file, force=True, disable_infer_genes=True,
        disable_infer_transcripts=True
):
    db = gffutils.create_db(
        annotation_file, dbfn=db_file, force=force, # change to True if we need to create a new db
        keep_order=True, merge_strategy='merge', sort_attribute_values=True,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts
    )
    return db


def main():
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "--annotation_file",
        dest="annotation_file",
        help="gff or gtf file",
        required=True
    )
    parser.add_argument(
        "--db_file",
        dest="db_file",
        help="database file output",
        required=True
    )
    parser.add_argument(
        "--force",
        action='store_true',
        help="if true, override any existing database file if it exists "
             "in the same path as specified output.",
        required=False,
        default=False
    )
    parser.add_argument(
        "--disable_infer_genes",
        action='store_true',
        help="If True, gffutils won't infer genes (default: False).",
        required=False,
        default=False
    )
    parser.add_argument(
        "--disable_infer_transcripts",
        action='store_true',
        help="If True, gffutils won't infer transcripts (default: False).",
        required=False,
        default=False
    )
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    annotation_file = args.annotation_file
    db_file = args.db_file
    force = args.force
    disable_infer_genes = args.disable_infer_genes
    disable_infer_transcripts = args.disable_infer_transcripts

    build_db(annotation_file, db_file, force, disable_infer_genes, disable_infer_transcripts)
if __name__ == "__main__":
    main()