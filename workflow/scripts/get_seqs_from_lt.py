import os
import pandas as pd
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
import argparse
import glob
import re


def main():
    locustags, lt, cds_faa, outfile = parsing_arguments()
    ## if locustags file exists, then only give locustags to filter_cds_faa
    if locustags != None:
        locustags = os.path.abspath(locustags)
        filtered_seqs, locus_tags = filter_cds_faa(cds_faa, locus_tags = locustags)
        ## if I want to recup the locus_tags for which I didn't have a match:
        descriptions = [seq.description for seq in filtered_seqs.values()]
        descriptions_string = " ".join(descriptions)
        missing_lts = [locus_tag for locus_tag in locus_tags if locus_tag not in descriptions_string]
        print(f"Missing locus_tags ({len(missing_lts)}) in the fasta file: {missing_lts}")
        ## write filtered_seqs (which will be a SeqRecord)
        with open(f"{outfile}", mode='w') as out:
            SeqIO.write(filtered_seqs.values(), out, "fasta")
        print(f"There were {len(locus_tags)} locus_tags in the input file, and {len(filtered_seqs)} sequences were found in {cds_faa}, and written to {outfile}.")
    elif lt != None:
        filtered_seqs, locus_tags = filter_cds_faa(cds_faa, single_lt = lt)
        with open(f"{outfile}", mode='a') as out:
            SeqIO.write(filtered_seqs.values(), out, "fasta")
        print(f"There were {len(locus_tags)} locus_tags in the input file, and {len(filtered_seqs)} sequences were found in {cds_faa}, and written to {outfile}.")


def parsing_arguments():
    parser = argparse.ArgumentParser(
    prog='get_seqs_from_lt.py',
    description="""
                Extract sequences from fasta file containing CDSs, based on locus_tags. 
                If a file is given with several locus_tags (1 per line), will create one output file. 
                If used in a loop (with for example 1 locus_tag per CDS file), will append to the same output file.
                """)
    parser.add_argument('--locus-tags',
                        help="File containing one locus_tag to retrieve per line.", required = False, default = None)
    parser.add_argument('--single-lt', 
                         help = "Single locus_tag to retrieve from cds_faa file.", required = False, default = None)
    parser.add_argument('--cds',
                        help="Fasta file containing CDSs to extract (AA or DNA, doesn't matter).")
    parser.add_argument('--out', help = 'Name of output file containing target sequences.')
    args = vars(parser.parse_args())
    locustags = args["locus_tags"]
    lt = args["single_lt"]
    cds_faa = os.path.abspath(args["cds"])
    outfile = os.path.abspath(args["out"])
    return locustags, lt, cds_faa, outfile

def filter_cds_faa(cds_faa, locus_tags = None, single_lt = None):
    ## read locus_tags file and store in list (array) (filter removes empty lines, and splitlines removes the '\n' at the end of each line)
    ## the list() assures the returned object is a list and not a filter object.
    print(locus_tags)
    if locus_tags != None:
        tags = list(filter(None, open(locus_tags, mode = "r").read().splitlines()))
    elif single_lt != None:
        tags = [single_lt]
    else:
        print("No locus-tags or single locus-tag provided.")
    ## make fasta dictionary from cds.faa file
    targets = {}
    with open(cds_faa, mode='rt') as f:
        faadic = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        for key, value in faadic.items():
            ## return locus_tag (PIINxxxx) from tags if locus_tag is present in string 'value.description' from the cds_faa dictionary 
            #res = [locus_tag for locus_tag in tags if (locus_tag in value.description)]
            #res = [locus_tag for locus_tag in tags if (f"{locus_tag}]" in value.description)]
            res = [locus_tag for locus_tag in tags if (f"{locus_tag}" in key)]
            ## if the result (res) is not en empty list ([]), then add the key and value to the targets dictionary
            if res != []:
                targets[key] = value
    ## return the final targets dictionary
    return targets, tags

if __name__ == '__main__':
    main()
