from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils
from pathlib import Path
import time
import gzip
import re
import argparse
import glob
import os

def main():
    gff_dir = parsing_arguments()
    
    ## Make and write GFF dictionaries
    gff_dir_to_dict(gff_dir)

def parsing_arguments():
    parser = argparse.ArgumentParser(
    prog='make_gff_db.py',
    description='Make database files for GFF files, to make them faster to process later.')
    parser.add_argument('--gff-dir', 
                        help = "Path to the directory with the Species.gff(.gz) files (genome annotations). Will create db in same directory")

    args = vars(parser.parse_args())
    gff_dir = os.path.abspath(args["gff_dir"])

    return gff_dir

def gff_dir_to_dict(gff_dir):
    """
    This function takes a directory containing the strain.gff and/or strain.gff.gz, and writes database files in the same directory.
    If the database files already exist, they will be overwritten.
    """
    print(f"Creating GFF databases in {gff_dir}. If dbs already exist, will overwrite them. You have 30 seconds to cancel if you want.")
    #time.sleep(30)
    ## can't glob.glob(*.gff.*) because if db has already been created, then it will also match *.gff.db, and that will raise error
    gffs = glob.glob(f"{gff_dir}/*.gff.gz") + \
        glob.glob(f"{gff_dir}/*.gff")
    for gff in gffs:
        path = gff.split(".")[0]
        key = gff.split("/")[-1].split(".")[0] # strain number
        db = gffutils.create_db(data=gff,
                                #dbfn=":memory:",
                                dbfn=f"{path}.gff.db",
                                merge_strategy='create_unique',
                                disable_infer_transcripts=True,
                                disable_infer_genes=True,
                                force = True,
                                keep_order = True)
        print(f"Created GFF database for {key} in {path}.gff.db")

if __name__ == '__main__':
    main()