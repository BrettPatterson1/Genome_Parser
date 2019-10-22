import argparse
import csv
from Bio import Align
import sys


aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0


def score(seqA, seqB):
    return aligner.score(seqA, seqB)

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

    return parser.parse_args()


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
    data[0].append(header_append)
    for i in range(1, len(data)):
        data[i].append(score(sequence, data[i][num_to_score]))


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
