from argparse import ArgumentParser
import Annotator
import sys

GENE_PRIORITY = [
    ['protein_coding','CDS'],
    ['protein_coding','start_codon'],
    ['protein_coding','stop_codon'],
    ['protein_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['protein_coding','5UTR'],
    ['protein_coding','3UTR'],
    ['protein_coding','UNCLASSIFIED_UTR'],
    ['protein_coding','intron'],
    ['protein_coding','Selenocysteine'],
    ['non_coding','exon'],
    ['non_coding','intron'],
    ['non_coding','transcript'],
    ['non_coding','gene'],
    ['non_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['non_coding','3UTR'],
    ['non_coding','5UTR'],
    ['non_coding','UNCLASSIFIED_UTR'],
    ['non_coding','Selenocysteine'],
    ['non_coding','CDS'],  # shouldn't occur?
    ['non_coding','start_codon'],  # shouldn't occur?
    ['non_coding','stop_codon'],  # shouldn't occur?
    ['protein_coding', 'exon'],  # shouldn't occur?
    ['protein_coding', 'transcript'],  # shouldn't occur?
    ['protein_coding', 'gene'],  # shouldn't occur?
]

TRANSCRIPT_PRIORITY = [
    ['protein_coding','CDS'],
    ['protein_coding','start_codon'],
    ['protein_coding','stop_codon'],
    ['protein_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['protein_coding','5UTR'],
    ['protein_coding','3UTR'],
    ['protein_coding', 'UNCLASSIFIED_UTR'],
    ['protein_coding','intron'],
    ['protein_coding','Selenocysteine'],
    ['non_coding', 'exon'],
    ['non_coding', 'intron'],
    ['non_coding','transcript'],
    ['non_coding','gene'],
    ['non_coding','THREE_AND_FIVE_PRIME_UTR'],
    ['non_coding','3UTR'],
    ['non_coding','5UTR'],
    ['non_coding', 'UNCLASSIFIED_UTR'],
    ['non_coding','Selenocysteine'],
    ['non_coding','CDS'], # shouldn't occur?
    ['non_coding','start_codon'],  # shouldn't occur?
    ['non_coding','stop_codon'],  # shouldn't occur?
    ['protein_coding','exon'], # shouldn't occur?
    ['protein_coding','transcript'], # shouldn't occur?
    ['protein_coding','gene'], # shouldn't occur
]

GENE_ID = 'gene_id'
TRANSCRIPT_ID = 'transcript_id'


def parse(priority_file, delim=','):
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
        help="gtf database file create from gffutils",
        required=True
    )
    parser.add_argument(
        "--limit-chroms-to",
        dest="limit_chroms_to",
        required=False,
        nargs='+',
        default=[]
    )
    parser.add_argument(
        "--transcript-priority-file",
        dest="transcript_priority_file",
        required=False,
        default=None
    )
    parser.add_argument(
        "--gene-priority-file",
        dest="gene_priority_file",
        required=False,
        default=None
    )
    parser.add_argument(
        "--unstranded",
        dest="unstranded",
        required=False,
        action='store_true',
        default=False
    )
    parser.add_argument(
        "--append_chr",
        dest="append_chr",
        required=False,
        action='store_true',
        default=False
    )
    parser.add_argument(
        "--species",
        dest="species",
        required=False,
        default='hg19'
    )
    parser.add_argument(
        "--fuzzy",
        dest="fuzzy",
        required=False,
        type=int,
        default=0
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
    fuzzy = args.fuzzy

    if unstranded:
        stranded=False
    else:
        stranded=True

    if args.transcript_priority_file is not None:
        transcript_priority = parse(args.transcript_priority_file)
    else:
        transcript_priority = TRANSCRIPT_PRIORITY

    if args.gene_priority_file is not None:
        gene_priority = parse(args.gene_priority_file)
    else:
        gene_priority = GENE_PRIORITY


    Annotator.annotate(
        gtfdb_file, input_bed_file, output_annotated_file, stranded, chroms,
        transcript_priority, gene_priority, species, append_chr, fuzzy
    )

if __name__ == "__main__":
    main()