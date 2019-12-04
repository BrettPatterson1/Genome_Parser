import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import re
import csv


def write_to_file(forward, reverse):
    try:
        with open("translated.txt", "w") as fileMake:
            fileMake.write(forward)
            fileMake.write("\n\n\n\n\n\n\n\n\n\n\n\n\n")
            fileMake.write(reverse)
            fileMake.close()
    except FileNotFoundError as fnf_error:
        print(fnf_error)


#Assumes Fasta File format
def read_stdin():
    handle = sys.stdin
    full_sequence = []
    for record in SeqIO.parse(handle, "fasta"):
        full_sequence.append(str(record.seq))
    return "".join(full_sequence)

def translate(dna):
    forward_sequence = Seq(dna, IUPAC.unambiguous_dna)

    # Gets reverse complement
    reverse_compliment_sequence = forward_sequence.reverse_complement()

    return str(forward_sequence.translate()), str(reverse_compliment_sequence.translate())

def main():
    dna_as_string = read_stdin()
    forward, reverse = translate(dna_as_string)
    write_to_file(forward, reverse)


if __name__ == '__main__':
    main()