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

def main():
    outdir, genus, level, overwrite, db, repo, max_g = parsing_arguments()
    ## Take care of directories: make the directory, then go to it
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)
    info = Path("./assembly.info")
    info.unlink(missing_ok=True)
    ## Set the timestamp
    month = date.today().strftime("%b%Y")
    ## Download the assembly_summary file (if doesn't exist)
    if not os.path.exists(f"assembly_summary_{month}.txt"):
        print(f"###### Downloading latest assembly_summary.txt (bacterial) from NCBI FTP server as assembly_summary_{month}.txt")
        os.system(f"rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/{db}/{repo}/assembly_summary.txt ./assembly_summary_{db}_{month}.txt")
    else:
        print(f"##### assembly_summary (bacterial) for {month} already exists.")
    ## Make the dry_run.txt using ncbi_genome_download (ngd). Using ngd to download the genomes didn't work for me (kept getting disconnected).
    dry_run_tab = make_dry_run(date = month, repo = repo, genus = genus, level = level)
    ## Make ftp_list by looking for the GCF from the dry_run_tab, in assembly_summary
    print(f"Reading assembly_summary, and writing ftp paths for selected {genus} and {level}")
    ftp_list = filter_assembly(file = f"assembly_summary_{db}_{month}.txt", dry_run_tab = dry_run_tab, date = month)
    print(len(ftp_list), "genomes found for", genus)
    if len(ftp_list) >= 300 and max_g == math.inf:
        print(f"##### WARNING: More than 300 genomes found for {genus} ({len(ftp_list)}). Please check the assembly_summary file.")
        n_genomes = input("How many genomes do you want to download? (press ENTER to accept default of 300):" ) or "300"
        n_genomes = int(n_genomes)
    else:
        n_genomes = min(len(ftp_list), max_g)
    download_gcf(ftp_list, overwrite, n_genomes)
    make_info_file(f"assembly_summary_{db}_{month}.txt", ftp_list, n_genomes)

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='download_genomes.py',
        description='Downloads GCF records from the ncbi ftp server for bacteria')

    parser.add_argument('--outdir', default = "./genomes", help = "Path to output directory. Defaults to current_dir/genomes.")
    parser.add_argument('--genus', help = "Specify genus of interest", required=True)
    parser.add_argument('--level-of-assembly', dest='level', help = "Level of assembly for download", choices = ['complete', 'chromosome', 'contig', 'scaffold', 'all'], default = "Complete Genome")
    parser.add_argument('--overwrite', help = "Overwrite assembly files if they already exist", default = "True", choices = ["True", "False"])
    parser.add_argument('--database', help = "Download either from 'refseq' or 'genbank'", default = 'refseq', choices = ['refseq', 'genbank'])
    parser.add_argument('--repository', help = 'Which organisms to download', choices = ['bacteria', 'fungi', 'invertebrate', 'metagenomes', 'other', 'plant', 
                                                                                      'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral'])
    parser.add_argument('--max-genomes', dest = "max", help = "Maximum number of genomes to download", default=math.inf, type = int)
    args = vars(parser.parse_args())
    outdir = args["outdir"]
    genus = args["genus"]
    level = args["level"]
    db = args["database"]
    repo = args["repository"]
    overwrite = args["overwrite"]
    max_g = args["max"]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} did not exist, creating it.")
    return outdir, genus, level, overwrite, db, repo, max_g

def make_dry_run(date, repo, genus, level):
    ## if dry_run.txt already exists for this date, use it
    if not os.path.exists(f"dry_run_{date}.txt"):
        print(f"dry_run_{date}.txt doesn't exist, making it.")
        ## Make the dry_run.txt using ncbi_genome_download (ngd). Using ngd to download the genomes didn't work for me (kept getting disconnected).
        cmd = f"ncbi-genome-download {repo} --genera {genus} --assembly-levels {level} --refseq-categories reference --output-folder ./ --dry-run > ./dry_run_{date}.txt"
        subprocess.run(cmd, shell = True) # tab-separated output, each genome on one line
    else:
        print(f"dry_run_{date}.txt already exists, using it.")
    with open(f"dry_run_{date}.txt", "r") as f:
        dry_run_tab = pd.read_csv(f, sep = "\t", header = None, skiprows = 1) # skip first line
    dry_run_tab.columns = ["gcf", "species", "strain"]
    return(dry_run_tab)

def filter_assembly(file, dry_run_tab, date):
    ## Load assembly summary file (assembly_file)
    assembly_file = pd.read_csv(file, sep="\t", header=1, low_memory=False)
    pd.options.display.max_colwidth = 500
    ## Filter assembly_file by looking for the GCF from dry_run_tab
    assembly_file_filt = assembly_file[assembly_file["#assembly_accession"].isin(dry_run_tab["gcf"])]
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
