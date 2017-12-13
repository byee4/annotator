#!/usr/bin/env python

# transitionning to python2/python3 support

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
from annotator import annotate_bed
import sys

GENE_PRIORITY = [
    ['protein_coding','CDS'],
    ['protein_coding','start_codon'],
    ['protein_coding','stop_codon'],
    ['protein_coding','5utr_and_3utr'],
    ['protein_coding','5utr'],
    ['protein_coding','3utr'],
    ['protein_coding','unclassified_utr'],
    ['protein_coding','intron'],
    ['protein_coding','Selenocysteine'],
    ['non_coding','exon'],
    ['non_coding','intron'],
    ['non_coding','transcript'],
    ['non_coding','gene'],
    ['non_coding','Selenocysteine'],
]

TRANSCRIPT_PRIORITY = [
    ['protein_coding','CDS'],
    ['protein_coding','start_codon'],
    ['protein_coding','stop_codon'],
    ['protein_coding','5utr_and_3utr'],
    ['protein_coding','5utr'],
    ['protein_coding','3utr'],
    ['protein_coding', 'unclassified_utr'],
    ['protein_coding','intron'],
    ['protein_coding','Selenocysteine'],
    ['non_coding', 'exon'],
    ['non_coding', 'intron'],
    ['non_coding','transcript'],
    ['non_coding','gene'],
    ['non_coding','Selenocysteine'],
]

def parse_annotation_priority(priority_file, delim=','):
    """
    From a comma delimited file, return a list of lists
    describing the priority with which to annotate stuff.

    :param priority_file:
    :return:
    """
    priority = []
    with open(priority_file, 'r') as f:
        for line in f:
            line = line.rstrip().split(delim)
            priority.append([line[0], line[1]])
    return priority

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
        help="input bed6 file",
        required=True
    )
    parser.add_argument(
        "--gtfdb",
        dest="gtfdb",
        help="gtf database file create from gffutils (Currently only "
             "GTF files are allowed)",
        required=True
    )
    parser.add_argument(
        "--transcript-priority-file",
        dest="transcript_priority_file",
        required=False,
        help='For each transcript (grouped by genes), return the transcript'
             ' using this priority',
        default=None
    )
    parser.add_argument(
        "--gene-priority-file",
        dest="gene_priority_file",
        required=False,
        help='For each gene, return the gene region using this priority',
        default=None
    )
    parser.add_argument(
        "--unstranded",
        dest="unstranded",
        required=False,
        action='store_true',
        help='ALLOW unstranded - will still search for correctly'
             ' stranded features first (searches positive first if a BED3)',
        default=False
    )
    parser.add_argument(
        "--append-chr",
        dest="append_chr",
        required=False,
        action='store_true',
        help='If your bedfile is formatted as something like [chrI] but the '
             'gtf database expects [I], we can retroactively append [chr] '
             'for you.',
        default=False
    )
    parser.add_argument(
        "--species",
        dest="species",
        required=False,
        help='sets default GTF nomenclature to either [ce11] or [hg19/mm10]',
        default='hg19'
    )
    parser.add_argument(
        "--limit-chroms-to",
        dest="limit_chroms_to",
        required=False,
        help="limit to annotating to these chromosomes "
             "only, saves time/memory",
        nargs='+',
        default=[]
    )
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    input_bed_file = args.input
    output_annotated_file = args.output
    gtfdb_file = args.gtfdb
    chroms = args.limit_chroms_to
    unstranded = args.unstranded
    species = args.species
    append_chr = args.append_chr
    cores = 1  # TODO: implement later - unnecessary now
    fuzzy = 0  # TODO: implement later - unnecessary now

    ### Flipping this terminology to make it easier to understand
    if unstranded:
        stranded=False
    else:
        stranded=True

    ### If appropriate priorities are unassigned, assign to defaults
    if args.transcript_priority_file is not None:
        transcript_priority = parse_annotation_priority(args.transcript_priority_file)
    else:
        transcript_priority = TRANSCRIPT_PRIORITY

    if args.gene_priority_file is not None:
        gene_priority = parse_annotation_priority(args.gene_priority_file)
    else:
        gene_priority = GENE_PRIORITY

    ### Call main function
    annotate_bed.annotate_bed(
        gtfdb_file, input_bed_file, output_annotated_file, stranded, chroms,
        transcript_priority, gene_priority, species, append_chr, fuzzy, cores
    )

if __name__ == "__main__":
    main()