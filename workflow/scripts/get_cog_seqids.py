import pandas as pd
import os 
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import glob
import gzip

def main():
    cog_dir, outfile = parsing_arguments()
    with open(outfile, mode = 'w') as out:
        out.write("cog;id\n")
    os.chdir(cog_dir)
    fastaList = glob.glob(f"*.fasta")
    for fasta in fastaList:
        cog = fasta.replace(".fasta", "")
        with open(fasta, mode = 'rt') as f:
            fastdic = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        for key in fastdic.keys():
            with open(outfile, mode = 'a') as out:
                out.write(f"{cog};{key}\n")

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='get_cog_seqids.py',
        description='Analyses COGs and outputs dataframe containing their names and the included seqids')

    parser.add_argument('--path', 
                        help="Path to the directory containing the COG files.")
    parser.add_argument('--outfile',
                        help="Name of the outfile to be created, a csv dataframe, separated by ';' containing the COG name for each seqid.")

    args = vars(parser.parse_args())
    cog_dir = os.path.abspath(args["path"])
    outfile = os.path.abspath(args["outfile"])

    return cog_dir, outfile


if __name__ == '__main__':
    main()