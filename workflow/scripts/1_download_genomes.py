import os
import argparse
from pathlib import Path
import gzip
import glob
from datetime import date
import pandas as pd
import re
import math
#import ncbi_genome_download as ngd
import subprocess
import random

def main():
    outdir, genus, level, overwrite, db, repo, max_g, req = parsing_arguments()
    print(f"Required genomes: {req}")
    ## Update max_g if req is not empty (i.e. if there are required genomes)
    if len(req) > 0 and max_g != math.inf:
        print(f"Counted {len(req)} required genomes, and max genomes limit set to {max_g}. {len(req)}/{max_g} will be your required genomes, and the rest will be 'randomly;' selected.")
        if len(req) > max_g:
            print(f"##### WARNING: You have more required genomes ({len(req)}) than the max-genomes limit ({max_g + len(req)}). Setting max-genomes to {len(req)}.")
            max_g = len(req)
        else:
            print(f"After accounting for required genomes, you can download {max_g} additional genomes.")
    
    ## Replace db with list of "genbank","refseq" if it is set to "both"
    if db == "both":
        dbs = ["genbank", "refseq"]
    else:
        dbs = [db]

    ## Take care of directories: make the directory, then go to it
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)
    info = Path("./assembly.info")
    info.unlink(missing_ok=True)
    
    ## Set the timestamp
    month = date.today().strftime("%b%Y")
    
    ## Download the assembly_summary file (if doesn't exist)
    # takes care of "genbank", "refseq" or "both" now. If both, concatenates the individual files
    download_assembly_summary(databases = dbs, repository = repo, date = month)

    ## Make the dry_run.txt using ncbi_genome_download (ngd). Using ngd to download the genomes didn't work for me (kept getting disconnected).
    # takes care of "genbank", "refseq" or "both" now. If both, concatenates the individual files
    # I'm still using ngd, because it enables me to use only "reference" genomes, de-duplicating the data quite a bit.
    dry_run_tab = make_dry_run(date = month, repo = repo, genus = genus, level = level, databases = dbs)
    
    ## Make ftp_list by looking for the GCF/GCA from the dry_run_tab and required genomes, in assembly_summary
    print(f"Reading assembly_summary, and writing ftp paths for selected {genus} and {level}")
    ftp_list = filter_assembly(file = f"assembly_summary_{db}_{month}.txt", dry_run_tab = dry_run_tab, date = month, required_genomes = req) # if db= "both", this will work, as I didn't replace db, but created a dbs list
    print(len(ftp_list), "genomes found for", genus)

    ## Communicate how many of the required genomes were found in assembly_summary
    if len(req) > 0:
        print(f"Checking how many of the {len(req)} required genomes were found in assembly_summary_{db}_{month}.txt")
        found_req = [req_acc for req_acc in req if any(req_acc in ftp.split("/")[-1] for ftp in ftp_list)]
        print(f"Found {len(found_req)} out of the {len(req)} required genomes in assembly_summary_{db}_{month}.txt")
    
    ## Determine how many genomes to download, if more than 300 are found
    if len(ftp_list) >= 300 and max_g == math.inf:
        print(f"##### WARNING: More than 300 genomes found for {genus} ({len(ftp_list)}). Please check the assembly_summary file.")
        n_genomes = input("How many genomes do you want to download? (press ENTER to accept default of 300):" ) or "300"
        n_genomes = int(n_genomes)
    else:
        n_genomes = min(len(ftp_list), max_g)
    print(f"Downloading {n_genomes} genomes.")

    ## Shuffle the ftp_list
    # First extract required genomes as req_ftp list, and remove them from the ftp_list
    req_ftps = [ftp for ftp in ftp_list if any(req_acc in ftp.split("/")[-1] for req_acc in found_req)]
    shuffled_ftps = [ftp for ftp in ftp_list if ftp not in req_ftps]
    
    ## Shuffle the shuffled_ftps
    random.seed(42) # for reproducibility
    random.shuffle(shuffled_ftps)
    ## Add the required genomes again, at the start of the list
    shuffled_ftps = req_ftps + shuffled_ftps

    ## Download a random subset of genomes (which HAS TO INCLUDE THE REQUIRED GENOMES, IF ANY)
    download_gcf(ftp_list = shuffled_ftps, overwrite = overwrite, n_genomes = n_genomes)
    make_info_file(f"assembly_summary_{db}_{month}.txt", shuffled_ftps, n_genomes)

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='download_genomes.py',
        description='Downloads GCF records from the ncbi ftp server for bacteria')

    parser.add_argument('--outdir', default = "./genomes", help = "Path to output directory. Defaults to current_dir/genomes.")
    parser.add_argument('--genus', help = "Specify genus of interest", required=True)
    parser.add_argument('--level-of-assembly', dest='level', help = "Level of assembly for download", choices = ['complete', 'chromosome', 'contig', 'scaffold', 'all'], default = "complete")
    parser.add_argument('--overwrite', help = "Overwrite assembly files if they already exist", default = "True", choices = ["True", "False"])
    parser.add_argument('--database', help = "Download either from 'refseq' or 'genbank'", default = 'refseq', choices = ['refseq', 'genbank', 'both'])
    parser.add_argument('--repository', help = 'Which organisms to download', choices = ['bacteria', 'fungi', 'invertebrate', 'metagenomes', 'other', 'plant', 
                                                                                      'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral'], default = 'bacteria')
    parser.add_argument('--max-genomes', dest = "max", help = "Maximum number of genomes to download", default=math.inf, type = int)
    parser.add_argument('--required-genomes', help = "List of required genomes (GCF or GCA IDs) to download, even if they exceed the max-genomes limit", nargs = '*', default = [])
    args = vars(parser.parse_args())
    outdir = args["outdir"]
    genus = args["genus"]
    level = args["level"]
    db = args["database"]
    repo = args["repository"]
    overwrite = args["overwrite"]
    max_g = args["max"]
    req = args["required_genomes"]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} did not exist, creating it.")
    return outdir, genus, level, overwrite, db, repo, max_g, req

def download_assembly_summary(databases, repository, date):
    for db in databases:
        month = date
        if not os.path.exists(f"assembly_summary_{db}_{month}.txt"):
            print(f"###### Downloading latest assembly_summary.txt (for {db} - {repository}) from NCBI FTP server as assembly_summary_{db}_{month}.txt")
            os.system(f"rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/{db}/{repository}/assembly_summary.txt ./assembly_summary_{db}_{month}.txt")
        else:
            print(f"##### assembly_summary (for {db} - {repository}) for {month} already exists.")
    # If 1 summary for refseq and one for genbank, concatenate them
    if len(databases) == 2:
        os.system(f"( cat assembly_summary_refseq_{month}.txt \"\n\"; sed '1,2d' assembly_summary_genbank_{month}.txt; echo \"\n\";) > assembly_summary_both_{month}.txt")
    return

def make_dry_run(date, repo, databases, genus, level):
    for db in databases:
        ## if dry_run.txt already exists for this date, use it
        if not os.path.exists(f"dry_run_{db}_{date}.txt"):
            print(f"dry_run_{db}_{date}.txt doesn't exist, making it.")
            ## Make the dry_run.txt using ncbi_genome_download (ngd). Using ngd to download the genomes didn't work for me (kept getting disconnected).
            cmd = f"ncbi-genome-download {repo} --genera {genus} --section {db} --assembly-levels {level} --refseq-categories reference --output-folder ./ --dry-run > ./dry_run_{db}_{date}.txt"
            print("Running command:", cmd)
            subprocess.run(cmd, shell = True) # tab-separated output, each genome on one line
        else:
            print(f"dry_run_{db}_{date}.txt already exists, using it.")
    ## Concatenate both, if len(databases) == 2
    if len(databases) == 2:
        os.system(f"( cat dry_run_refseq_{date}.txt \"\n\"; sed '1d' dry_run_genbank_{date}.txt; echo \"\n\";) > dry_run_both_{date}.txt")
        db = "both"
    else:
        db = databases[0]
    with open(f"dry_run_{db}_{date}.txt", "r") as f:
        dry_run_tab = pd.read_csv(f, sep = "\t", header = None, skiprows = 1)
    dry_run_tab.columns = ["acc", "species", "strain"]
    return(dry_run_tab)

def filter_assembly(file, dry_run_tab, date, required_genomes):
    ## Load assembly summary file (assembly_file)
    assembly_file = pd.read_csv(file, sep="\t", header=1, low_memory=False)
    pd.options.display.max_colwidth = 500
    ## Get all accessions from dry_run_tab AND required_genomes
    accessions = dry_run_tab["acc"].tolist() + required_genomes
    accessions = list(set(accessions)) # remove duplicates, if any
    ## Filter assembly_file by looking for the GCF from dry_run_tab
    assembly_file_filt = assembly_file[assembly_file["#assembly_accession"].isin(accessions)] # CHECK THAT IT WORKS
    ## Take only the ftp_path column, and modify as needed
    ftp_path = assembly_file_filt["ftp_path"].to_string(index=False)
    ftp_path = re.sub("https://", "", ftp_path)
    ftp_list = list(map(str.strip, ftp_path.split("\n")))
    with open(f"ftpdirpaths_{date}.txt", mode ='wt', encoding = "utf-8") as myftp:
        myftp.write('\n'.join(ftp_list))
    with open(f"assembly_summary_filtered_{date}.tsv", mode = 'wt', encoding = "utf-8") as myasm:
        assembly_file_filt.to_csv(myasm, sep = "\t", index = False)
    return ftp_list

def download_gcf(ftp_list, overwrite, n_genomes):
    ## Loop over range(n_genomes), to download the first n_genomes from shuffled_ftps
    for i in range(n_genomes):
        gcf = ftp_list[i].split("/")[-1]
        if overwrite == "True":
            print(f"Downloading assembly {gcf}")
            os.system(
                f"rsync --copy-links --times --verbose -r --keep-dirlinks rsync://{ftp_list[i]} ./")
        else:
            if os.path.exists(gcf):
                print(
                    f"Directory for assembly {gcf} already exists, skipping...")
                continue
            else:
                print(f"Downloading assembly {gcf}")
                os.system(
                    f"rsync --copy-links --times --verbose -r --keep-dirlinks rsync://{ftp_list[i]} ./")

def make_info_file(file, ftp_list, n_genomes):
    #for i in range(10):
    for i in range(n_genomes):
        assembly_file = pd.read_csv(file, sep="\t", header=1, low_memory=False)
        gcf = ftp_list[i].split("/")[-1]
        print(gcf)
        gbff = glob.glob(f"./{gcf}/*gbff*")[0]
        organism = assembly_file[assembly_file["ftp_path"].str.contains(gcf)]["organism_name"].to_string(index = False)
        if gbff.endswith(".gz"):
            with gzip.open(gbff, mode = 'r') as zipf:
                lt, olt = get_lt_olt(zipf, zip = True)
        else:
            with open(gbff, mode = 'r') as f:
                lt, olt = get_lt_olt(f, zip = False)
        with open("assembly.info", 'a') as info:
            info.write(f"{i}\t{gcf}\t{organism}\t{lt}\t{olt}\n")

def get_lt_olt(filehandle, zip = True):
    if zip == True:
        content = filehandle.read().decode("utf-8")
    else:
        content = filehandle.read()
    locustag = re.search(
        "(?i)/locus_tag=\"([A-Z0-9]*)_*", content)
    if locustag is not None:
        lt = locustag.group(1)
    else:
        lt = ""
    old_locustag = re.search(
        "(?i)/old_locus_tag=\"([A-Z0-9]*)_*", content)
    if old_locustag is not None:
        olt = old_locustag.group(1)
    else:
        olt = ""
    return(lt, olt)

if __name__ == '__main__':
    main()
