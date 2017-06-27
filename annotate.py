from argparse import ArgumentParser
import Annotator
import sys

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
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    input_bed_file = args.input
    output_annotated_file = args.output
    gtfdb_file = args.gtfdb
    chroms = args.limit_chroms_to

    Annotator.annotate(
        gtfdb_file, input_bed_file, output_annotated_file, chroms
    )

if __name__ == "__main__":
    main()