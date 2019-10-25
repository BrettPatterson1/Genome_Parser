import argparse
import csv
from Bio import Align
from Bio import pairwise2
import sys


aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0


MATCH_SCORE = 2
MISMATCH_SCORE = -1
OPEN_GAP_SCORE = -0.5
EXTEND_GAP_SCORE = -0.1

def score(seqA, seqB):
    return aligner.score(seqA, seqB)

def get_score_and_sequences(seqA, seqB):
    alignments = aligner.align(seqA, seqB)
    maximum = max(alignments)
    as_string = str(maximum)
    split = as_string.split("\n")
    found = split[0]
    original_with_gaps = split[2]
    return [maximum.score, found, original_with_gaps]
    # alignment = pairwise2.align.globalms(seqA, seqB, MATCH_SCORE, MISMATCH_SCORE, OPEN_GAP_SCORE, EXTEND_GAP_SCORE, one_alignment_only = True)

def get_arguments():
    desc = "Sorts a csv by match score for particular arguments"

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("sequence", help="sequence to compare", type=str)

    parser.add_argument('-e', '--export',
                        help="file name to export (no extension)." +
                             "Defaults to printing to stdout",
                        type=str)

    parser.add_argument('-f', '--file',
                        help='input csv to parse, defaults to stdin', type=str)

    typ_seq = parser.add_mutually_exclusive_group()

    typ_seq.add_argument('-aa', '--amino_acids',
                         help="parse amino acids", action="store_true")

    typ_seq.add_argument('-up', '--upstream',
                         help="parse upstream", action="store_true")

    typ_seq.add_argument('-nuc', '--nucleotides',
                         help="parse nucleotide sequence", action="store_true")

    args, unknown = parser.parse_known_args()
    return args


def read_csv(filename):
    try:
        with open(filename, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            list = []
            for row in csv_reader:
                if len(row) > 0:
                    list.append(row)
            return list
    except FileNotFoundError as fnf_error:
        print(fnf_error)


def match(header_append, num_to_score, sequence, data):
    header_extension = [header_append, "Closest Match Found", "Original Sequence With Gaps"]
    data[0].extend(header_extension)
    for i in range(1, len(data)):
        data[i].extend(get_score_and_sequences(sequence, data[i][num_to_score]))
        # data[i].append(score(sequence, data[i][num_to_score]))


def get_matches(data, args):
    if args.amino_acids:
        match(header_append="AA Score", num_to_score=4, sequence=args.sequence, data=data)
    elif args.upstream:
        match(header_append="Upstream Score", num_to_score=5, sequence=args.sequence, data=data)
    else:
        match(header_append="Nucleotide Score", num_to_score=3, sequence=args.sequence, data=data)


def sort_data(data):
    data[1:len(data)] = sorted(data[1:len(data)], key= lambda x:x[6], reverse=True)


def write_to_file(output_filename, data):
    try:
        with open(output_filename, "w") as fileMake:
            w = csv.writer(fileMake)
            w.writerows(data)
            fileMake.close()
    except FileNotFoundError as fnf_error:
        print(fnf_error)


def read_in():
    lines = sys.stdin.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace('\n','')
        lines[i] = eval(lines[i])
    return lines


def main():
    args = get_arguments()
    if args.file:
        data = read_csv(args.file)
    else:
        data = read_in()
    get_matches(data, args)
    sort_data(data)
    if args.export:
        write_to_file(args.export, data)
    else:
        for row in data:
            print(row)



if __name__ == '__main__':
    main()
