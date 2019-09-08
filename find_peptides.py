import os

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import re
import csv


# "Main" of the program, checks to see if the file exists and if it does, parses the forward and reverse complement
def main():
    args = get_arguments()
    output_filename = create_output_filename(args)
    output_as_fasta = args.fasta
    all_peptide_information = extract_data(args)
    write_to_file(output_filename, output_as_fasta, all_peptide_information)


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
    start_codons = args.start
    output_as_fasta = args.fasta
    minimum_peptide_length = args.minlen
    maximum_peptide_length = args.maxlen
    all_peptide_information = []
    try:
        with open(args.file, "r") as handle:
            # Reads in the sequence
            genome_sequence = SeqIO.read(handle, "fasta").seq
            # Sets the alphabet to DNA
            forward_sequence = Seq(str(genome_sequence), IUPAC.unambiguous_dna)

            # Gets reverse complement
            reverse_compliment_sequence = forward_sequence.reverse_complement()
            all_peptide_information.extend(
                find_all_possible_proteins(forward_sequence, "+", minimum_peptide_length, maximum_peptide_length,
                                           start_codons,
                                           output_as_fasta))
            all_peptide_information.extend(
                find_all_possible_proteins(reverse_compliment_sequence, '-', minimum_peptide_length,
                                           maximum_peptide_length, start_codons,
                                           output_as_fasta))
    except FileNotFoundError as fnf_error:
        print(fnf_error)
    return all_peptide_information


def create_output_filename(args):
    # This just returns the file directory to read - the extension so we can make a new file based on it
    # If this causes errors in the future, take a look at os.splitext
    # https://docs.python.org/2/library/os.path.html
    # output_filename = ((fileDir_r[::-1]).split(".", 1)[1][::-1] + "_parsed")
    if args.export:
        output_filename = args.export
    else:
        output_filename = (os.path.splitext(args.file)[0] + "_parsed")
    if args.fasta:
        output_filename = output_filename + ".fasta"
    else:
        output_filename = output_filename + ".csv"
    return output_filename


# Organizes given data into a csv-format
def create_csv_sequence(index, mark, length, dna_as_string, protein_as_string):
    return [index, mark, length, dna_as_string, protein_as_string]


# Organizes given data into a fasta-format
def create_fasta_sequence(dna):
    return SeqRecord(dna)


# If gets slow, open the file to append and append instead of storing everything to memory


# Finds and prints all of the proteins given a sequence and a forward/reverse complement tag
def find_all_possible_proteins(genome_sequence, direction_indicator, minimum_peptide_length, maximum_peptide_length,
                               start_codons_as_string,
                               output_as_fasta):
    # splits the start argument
    start_codons = extract_start_codons(start_codons_as_string)

    # Finds the indices of all start codons- regardless of reading frame
    start_codon_locations = get_all_start_codon_locations(genome_sequence, start_codons)

    all_peptide_information = []
    for index in start_codon_locations:
        # Tests end case to avoid array out of bounds exception
        if window_exceeds_genome_array_bounds(genome_sequence, index, maximum_peptide_length):
            dna_window = genome_sequence[index:]
            translated_window = dna_window.translate(to_stop=True)
            if protein_meets_length_specifications(minimum_peptide_length, maximum_peptide_length, translated_window):
                dna_up_to_stop_codon = dna_window[0:(len(translated_window) + 1) * 3]
                # print("index %i, %s, Length: %i, DNA: %s, Protein: %s"
                # (index, direction indicator, len(translated_window), str(dna_up_to_stop_codon), str(translated_window)))
                if output_as_fasta:
                    all_peptide_information.append(create_fasta_sequence(dna_up_to_stop_codon))
                else:
                    all_peptide_information.append(
                        create_csv_sequence(index, direction_indicator, len(translated_window),
                                            str(dna_up_to_stop_codon),
                                            str(translated_window)))

        else:
            # Creates a subsequence of max length + 1
            dna_window = genome_sequence[index:(index + 3 * (maximum_peptide_length + 1))]
            # Translates the subsequence and cuts off just before a stop codon
            translated_window = dna_window.translate(to_stop=True)
            # Checks to see if the translated sequence is inbetween the inputted max and min and prints if true
            if protein_meets_length_specifications(minimum_peptide_length, maximum_peptide_length, translated_window):
                dna_up_to_stop_codon = dna_window[0:(len(translated_window) + 1) * 3]
                # print("index %i, %s, Length: %i, DNA: %s, Protein: %s"
                # (index, toMark, len(translated_window), str(dna_up_to_stop_codon), str(translated_window)))
                if output_as_fasta:
                    all_peptide_information.append(create_fasta_sequence(dna_up_to_stop_codon))
                else:
                    all_peptide_information.append(
                        create_csv_sequence(index, direction_indicator, len(translated_window),
                                            str(dna_up_to_stop_codon),
                                            str(translated_window)))
    return all_peptide_information


def protein_meets_length_specifications(minimum_peptide_length, maximum_peptide_length, translated_window):
    return minimum_peptide_length <= len(translated_window) <= maximum_peptide_length


def window_exceeds_genome_array_bounds(genome_sequence, index, maximum_peptide_length):
    return (3 * (maximum_peptide_length + 1)) + index > len(genome_sequence)


def get_all_start_codon_locations(sequence, start_codons):
    locations = []
    for i in range(len(start_codons)):
        locations.extend([m.start() for m in re.finditer(start_codons[i], str(sequence))])
    return locations


def extract_start_codons(start_codons_as_string):
    starters = start_codons_as_string.split(",")
    for i in range(len(starters)):
        starters[i] = starters[i].strip()
    return starters


if __name__ == '__main__':
    main()
