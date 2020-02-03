import os
import sys
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
    output_as_fasta = args.fasta

    if args.anchor_location_file is not None:
        anchor_locations = get_anchor_locations(args.anchor_location_file)
    else:
        anchor_locations = []
    dna_as_strings = get_dna_as_string(args.fasta_file)
    all_peptide_information = extract_data(args, anchor_locations, dna_as_strings)
    upstream_bases = args.upstream
    if args.export:
        output_filename = create_output_filename(args)
        write_to_file(output_filename, output_as_fasta, all_peptide_information, upstream_bases)
    else:
        if output_as_fasta:
            print(all_peptide_information)
        else:
            print(['Index', 'Direction', 'Protein Length', 'DNA Sequence', 'Protein Sequence', "Upstream {} Bases".format(upstream_bases)])
            for row in all_peptide_information:
                print(row)


def get_dna_as_string(fasta_file):
    if fasta_file is None:
        return read_stdin()
    try:
        with open(fasta_file, "rU") as handle:
            # Reads in the sequence
            full_sequence = []
            for record in SeqIO.parse(handle, "fasta"):
                full_sequence.append(str(record.seq))
            return "".join(full_sequence)
    except FileNotFoundError as fnf_error:
        print(fnf_error)


#Assumes Fasta File format
def read_stdin():
    handle = sys.stdin
    full_sequence = []
    for record in SeqIO.parse(handle, "fasta"):
        full_sequence.append(str(record.seq))
    return "".join(full_sequence)


def get_anchor_locations(filename):
    f = open(filename, "r")
    locations_as_string = ""
    if f.mode == 'r':
        locations_as_string = f.read()
    locations = []
    separate_locations = locations_as_string.split(" ")
    for location in separate_locations:
        if location == "":
            continue
        else:
            no_space_location = location.strip()
            if no_space_location != "":
                locations.append(int(no_space_location))
    return locations

def get_arguments():
    desc = "Finds all possible peptides given a range of amino acids."

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-e', '--export',
                        help="file name to export (no extension)." +
                             "Default prints to stdout",
                        type=str)

    parser.add_argument('-max', '--maxlen', default=45,
                        help='Maximum length of amino acids (45 default)', type=int)

    parser.add_argument('-u', "--upstream", default=120,
                        help="Number of nucleotide bases to record upstream", type=int)

    parser.add_argument('-min', '--minlen', default=15,
                        help='Minimum length of amino acids (15 default)', type=int)

    parser.add_argument('-s', '--start', default="ATG",
                        help='start codons to user (separated by commas)- ATG default', type=str)

    parser.add_argument('-f', '--fasta_file',
                        help='input fasta file to parse, else reads from stdin and assumes fasta_file format', type=str)

    parser.add_argument('-a', '--anchor_location_file', default=None,
                        help='text file with anchor indices based on nucleotide location', type=str)

    parser.add_argument('-r', '--radius', default=50,
                        help='Maximum radius tolerance from anchor locations', type=int)

    typ = parser.add_mutually_exclusive_group()

    typ.add_argument('-c', '--csv',
                     help="Export file as csv (default)", action="store_true")

    typ.add_argument('-fa', '--fasta',
                     help="Export file as fasta", action="store_true")

    args, unknown = parser.parse_known_args()
    return args


def write_to_file(output_filename, output_as_fasta, all_peptide_information, upstream_bases):
    try:
        with open(output_filename, "w") as fileMake:
            if output_as_fasta:
                SeqIO.write(all_peptide_information, fileMake, "fasta")
            else:
                w = csv.writer(fileMake)
                w.writerow(['Index', 'Direction', 'Protein Length', 'DNA Sequence', 'Protein Sequence', "Upstream {} Bases".format(upstream_bases)])
                w.writerows(all_peptide_information)
            fileMake.close()
    except FileNotFoundError as fnf_error:
        print(fnf_error)


def extract_data(args, anchor_locations, genome_sequence_as_str):
    radius = args.radius
    start_codons = args.start
    output_as_fasta = args.fasta
    minimum_peptide_length = args.minlen
    maximum_peptide_length = args.maxlen
    all_peptide_information = []
    upstream_length = args.upstream
    # Reads in the sequence
    # Sets the alphabet to DNA
    forward_sequence = Seq(genome_sequence_as_str, IUPAC.unambiguous_dna)

    # Gets reverse complement
    reverse_compliment_sequence = forward_sequence.reverse_complement()
    all_peptide_information.extend(
        find_all_possible_proteins(forward_sequence, "+", minimum_peptide_length, maximum_peptide_length,
                                   start_codons, output_as_fasta, radius, anchor_locations, upstream_length))
    all_peptide_information.extend(
        find_all_possible_proteins(reverse_compliment_sequence, '-', minimum_peptide_length,
                                   maximum_peptide_length, start_codons,
                                   output_as_fasta, radius, anchor_locations, upstream_length))
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
def create_csv_sequence(index, mark, length, dna_as_string, protein_as_string, upstream_bases):
    return [index, mark, length, dna_as_string, protein_as_string, upstream_bases]


# Organizes given data into a fasta-format
def create_fasta_sequence(dna):
    return SeqRecord(dna)


# If gets slow, open the file to append and append instead of storing everything to memory


# Finds and prints all of the proteins given a sequence and a forward/reverse complement tag
def find_all_possible_proteins(genome_sequence, direction_indicator, minimum_peptide_length, maximum_peptide_length,
                               start_codons_as_string, output_as_fasta, radius, anchor_locations, upstream_length):
    # splits the start argument
    start_codons = extract_start_codons(start_codons_as_string)

    # Finds the indices of all start codons- regardless of reading frame
    start_codon_locations = get_all_start_codon_locations(genome_sequence, start_codons, anchor_locations, radius, direction_indicator)
    length_of_sequence = len(genome_sequence)
    all_peptide_information = []
    for index in start_codon_locations:
        if index - upstream_length < 0:
            upstream_bases = str(genome_sequence[:index])
        else:
            upstream_bases = str(genome_sequence[index - upstream_length: index])

        # Tests end case to avoid array out of bounds exception
        if window_exceeds_genome_array_bounds(length_of_sequence, index, maximum_peptide_length):
            dna_window = genome_sequence[index:]
            len_window = len(dna_window)
            # Use this to make multiple of three
            dna_window = dna_window[:(len_window // 3) * 3]
            translated_window = dna_window.translate(to_stop=True)
            # Means never reached stop codon
            if len(translated_window) * 3 == len(dna_window):
                continue
            if protein_meets_length_specifications(minimum_peptide_length, maximum_peptide_length, translated_window):
                dna_up_to_stop_codon = dna_window[:(len(translated_window) + 1) * 3]
                # print("index %i, %s, Length: %i, DNA: %s, Protein: %s"
                # (index, direction indicator, len(translated_window), str(dna_up_to_stop_codon), str(translated_window)))
                if output_as_fasta:
                    all_peptide_information.append(create_fasta_sequence(dna_up_to_stop_codon))
                else:
                    corrected_index = index
                    if direction_indicator == "-":
                        corrected_index = length_of_sequence - index - 1
                    all_peptide_information.append(
                        create_csv_sequence(corrected_index, direction_indicator, len(translated_window),
                                            str(dna_up_to_stop_codon),
                                            str(translated_window), upstream_bases))

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
                    corrected_index = index
                    if direction_indicator == "-":
                        corrected_index = length_of_sequence - index - 1
                    all_peptide_information.append(
                        create_csv_sequence(corrected_index, direction_indicator, len(translated_window),
                                            str(dna_up_to_stop_codon),
                                            str(translated_window), upstream_bases))
    return all_peptide_information


def protein_meets_length_specifications(minimum_peptide_length, maximum_peptide_length, translated_window):
    return minimum_peptide_length <= len(translated_window) <= maximum_peptide_length


def window_exceeds_genome_array_bounds(len_genome_sequence, index, maximum_peptide_length):
    return (3 * (maximum_peptide_length + 1)) + index > len_genome_sequence


def get_all_start_codon_locations(sequence, start_codons, anchor_locations, radius, direction_indicator):
    all_locations = []
    for i in range(len(start_codons)):
        all_locations.extend([m.start() for m in re.finditer(start_codons[i], str(sequence))])
    if len(anchor_locations) == 0:
        return all_locations

    if direction_indicator == "-":
        length_of_sequence = len(sequence)
        for i in range(len(anchor_locations)):
            anchor_locations[i] = length_of_sequence - 1 - anchor_locations[i]
    filtered_locations = []
    for location in all_locations:
        for anchor in anchor_locations:
            if anchor - radius <= location <= anchor + radius:
                filtered_locations.append(location)
                break
    return filtered_locations


def extract_start_codons(start_codons_as_string):
    starters = start_codons_as_string.split(",")
    for i in range(len(starters)):
        starters[i] = starters[i].strip()
    return starters


if __name__ == '__main__':
    main()
