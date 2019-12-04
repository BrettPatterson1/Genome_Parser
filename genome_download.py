import argparse
from Bio import Entrez
import os
import urllib.request
import io
import gzip

Entrez.email = "brettsp2@illinois.edu"

def get_arguments():
    desc = "Downloads NCBI's Genome Files using Biopython Entrez"

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("species_file", help="file with names of species to download, one on each line", type=str)

    args, unknown = parser.parse_known_args()
    return args

def make_dir(path):
    try:
        os.mkdir(path)
    except OSError:
        print("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s " % path)


def download_species(species):
    handle = Entrez.esearch(db="assembly", term=(species + "[ORGN] AND latest[SB]"))
    record = Entrez.read(handle)
    handle.close()
    print(record)
    id_list = record['IdList']
    for seq_id in id_list:
        handle = Entrez.efetch(db="assembly", id=seq_id, rettype="fasta", retmode="text")
        seq = Entrez.read(handle)
        file = open('als.fasta', 'w')
        file.write(seq.rstrip("\n"))


def download_all_species(species):
    for name in species:
        saved_path = os.getcwd()
        path = os.path.join(saved_path, name)
        make_dir(path)
        os.chdir(path)
        get_assemblies(term=(name + "[ORGN] AND latest[SB]"), download=True)
        os.chdir(saved_path)


def main():
    args = get_arguments()
    species = get_species(args.species_file)
    download_all_species(species)


def get_species(filename):
    try:
        with open(filename, "r") as f:
            species = [line.rstrip('\n') for line in f]
            return species
    except FileNotFoundError as fnf_error:
        print(fnf_error)


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def get_assemblies(term, download=True):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = (url + "/" + label +'_genomic.fna.gz')
        print(link)
        links.append(link)
        if download == True:
            #download link
            filename = f'{label}.fna.gz'
            out_file_path = filename[:-3]
            response = urllib.request.urlopen(link)
            compressed_file = io.BytesIO(response.read())
            decompressed_file = gzip.GzipFile(fileobj=compressed_file)

            with open(out_file_path, "w") as outfile:
                outfile.write(decompressed_file.read().decode('utf-8'))
    return links


if __name__ == '__main__':
    main()