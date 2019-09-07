import os

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import re
import csv

isFasta = False
isCsv = False


def get_arguments():
    desc = "Finds all possible peptides given a range of amino acids."

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-e', '--export',
                        help="file name to export (no extension)." +
                             "Default is the file name with extension '_parsed'",
                        type=str)

    parser.add_argument('-max', '--maxlen', default=45,
                        help='Maximum length of amino acids (45 default)', type=int)

    parser.add_argument('-min', '--minlen', default=15,
                        help='Minimum length of amino acids (15 default)', type=int)

    parser.add_argument('-s', '--start', default="ATG",
                        help='start codons to user (separated by commas)- ATG default', type=str)

    parser.add_argument('file',
                        help='input file to parse', type=str)

    typ = parser.add_mutually_exclusive_group()

    typ.add_argument('-c', '--csv',
                     help="Export file as csv (default)", action="store_true")

    typ.add_argument('-f', '--fasta',
                     help="Export file as fasta", action="store_true")

    return parser.parse_args()

#"Main" of the program, checks to see if the file exists and if it does, parses the forward and reverse complement
def Main():
    args = get_arguments()
    output_filename = create_output_filename(args)
    output_as_fasta = args.fasta
    all_peptide_information = extract_data(args)
    write_to_file(output_filename, output_as_fasta, all_peptide_information)


def write_to_file(output_filename, output_as_fasta, all_peptide_information):
    try:
        with open(output_filename, "w") as fileMake:
            if output_as_fasta:
                SeqIO.write(all_peptide_information, fileMake, "fasta")
            else:
                w = csv.writer(fileMake)
                w.writerow(['Index', 'Direction', 'Protein Length', 'DNA Sequence', 'Protein Sequence'])
                w.writerows(all_peptide_information)
            fileMake.close()
    except FileNotFoundError as fnf_error:
        print(fnf_error)


def extract_data(args):
    minimum_peptide_length = args.minlen
    maximum_peptide_length = args.maxlen
    all_peptide_information = []
    try:
        with open(args.file, "r") as handle:
            # Reads in the sequence
            record_temp = SeqIO.read(handle, "fasta").seq
            # Sets the alphabet to DNA
            record = Seq(str(record_temp), IUPAC.unambiguous_dna)

            # Gets reverse complement
            record_rc = record.reverse_complement()
            all_peptide_information.extend(find(record, "+", minimum_peptide_length, maximum_peptide_length, args.start, args.fasta))
            all_peptide_information.extend(find(record_rc, '-', minimum_peptide_length, maximum_peptide_length, args.start, args.fasta))
    except FileNotFoundError as fnf_error:
        print(fnf_error)
    return all_peptide_information


def create_output_filename(args):
    # This just returns the file directory to read - the extension so we can make a new file based on it
    # If this causes errors in the future, take a look at os.splitext
    # https://docs.python.org/2/library/os.path.html
    # fileDir_w = ((fileDir_r[::-1]).split(".", 1)[1][::-1] + "_parsed")
    if args.export:
        fileDir_w = args.export
    else:
        fileDir_w = (os.path.splitext(args.file)[0] + "_parsed")
    if args.fasta:
        fileDir_w = fileDir_w + ".fasta"
    else:
        fileDir_w = fileDir_w + ".csv"
    return fileDir_w


# Organizes given data into a csv-format
def csvSequence(index, mark, length, DNAAsString, proteinAsString) :
    return [index, mark, length, DNAAsString, proteinAsString]


# Organizes given data into a fasta-format
def fastaSequence(DNA) :
    return SeqRecord(DNA)

# If gets slow, open the file to append and append instead of storing everything to memory


# Finds and prints all of the proteins given a sequence and a forward/reverse complement tag
def find(sequence, toMark, min, max, start, isFasta):

    # splits the start argument
    starters = start.split(",")
    for i in range(len(starters)):
        starters[i] = starters[i].strip()

    count = 0
    # Finds the indices of all start codons- regardless of reading frame
    locations = [m.start() for m in re.finditer(starters[0], str(sequence))]
    for i in range(1, len(starters)):
        locations.extend([m.start() for m in re.finditer(starters[i], str(sequence))])

    toRet = []
    for index in locations:
        # Tests end case to avoid array out of bounds exception
        if (len(sequence) - index < (3*(max + 1))):
            toTest = sequence[index:]
            toTest_tr = toTest.translate(to_stop=True)
            if (min <= len(toTest_tr) and len(toTest_tr) <= max):
                toTest_dna = toTest[0:len(toTest_tr)*3 + 3]
                # print("index %i, %s, Length: %i, DNA: %s, Protein: %s"
                # (index, toMark, len(toTest_tr), str(toTest_dna), str(toTest_tr)))
                count += 1
                if (isFasta):
                    toRet.append(fastaSequence(toTest_dna))
                else:
                    toRet.append(csvSequence(index, toMark, len(toTest_tr), str(toTest_dna), str(toTest_tr)))

        else:
            # Creates a subsequence of max length max+1
            toTest = sequence[index:(index + 3*(max + 1))]
            # Translates the subsequence and cuts off just before a stop codon
            toTest_tr = toTest.translate(to_stop=True)
            # Checks to see if the translated sequence is inbetween the inputted max and min and prints if true
            if (min <= len(toTest_tr) and len(toTest_tr) <= max):
                toTest_dna = toTest[0:len(toTest_tr)*3 + 3]
                # print("index %i, %s, Length: %i, DNA: %s, Protein: %s"
                # (index, toMark, len(toTest_tr), str(toTest_dna), str(toTest_tr)))
                count += 1
                if (isFasta):
                    toRet.append(fastaSequence(toTest_dna))
                else:
                    toRet.append(csvSequence(index, toMark, len(toTest_tr), str(toTest_dna), str(toTest_tr)))
    return toRet


if __name__ == '__main__':
    Main()

