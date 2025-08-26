from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re
import argparse
import os
import glob
import gzip

def main():
    genome_dir, infofile, logfile = parsing_arguments()
    with open(f"{logfile}", 'w') as log:
        log.write(f"""### This log file contains the countsof: 
            # CDSs (cds_from_genomic), 
            # Proteins (protein.faa), 
            # CDSs that are proteins (cds.faa created with this program),
            # Paralog pairs of proteins (CDSs that are proteins - Proteins)
            GCF\tCDS count\tProteins\tCDS-WP\tParalogPairs
            """)
    ## Change directory to genome_dir
    os.chdir(genome_dir)
    curr_wd = os.getcwd()

    ## Read the assembly.info file
    info = pd.read_csv(infofile, sep = "\t", header = None)
    ## Loop over each line (each strain/assembly information)
    for i in range(len(info)):
        strain, gcf = info.loc[i, 0:1]
        os.chdir(gcf)
        print(f"Working on {strain} ({gcf})")
        n_cds, n_prot, n_cds_prot, lt2number, number2wp = write_cds_faa(strain)
        with open(logfile, 'a') as log:
            log.write(f"{gcf}\t{n_cds}\t{n_prot}\t{n_cds_prot}\t{n_cds_prot - n_prot}\n")
        write_cds_gff(strain, lt2number, number2wp)
        write_correspondance_table(lt2number, number2wp, strain)

        os.chdir(curr_wd)

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='2_make_cds_faa_gff.py',
        description='Re-orders the protein.faa record for downloaded genomes, according to the order found in cds_from_genomic.fna')

    parser.add_argument('--genome_dir', default = "./genomes",
                        help="Path to the directory containing assembly directories (GCF_somethingsomething/).")
    parser.add_argument('--infofile', default = "./assembly.info",
                        help="Path to the assembly.info file that was created when downloading genomes.")
    parser.add_argument('--logfile', default = "./logfile", 
                        help="Log file to show if index.cds.faa has been created")

    args = vars(parser.parse_args())
    genome_dir = os.path.abspath(args["genome_dir"])
    infofile = os.path.abspath(args["infofile"])
    logfile = os.path.abspath(args["logfile"])

    return genome_dir, infofile, logfile

def fasta_gz_to_dict(pattern):
    """
    This function takes a filehandle of a fasta file, then reads the content and stores it in a dictionary, which it returns.
    It can handle gz zipped files, as well as unzipped fasta files, and doesn't care about the used alphabet.
    """
    fh = glob.glob(pattern)[0]
    if fh.endswith(".gz"):
        with gzip.open(fh, mode='rt') as zipf:
            content = SeqIO.to_dict(SeqIO.parse(zipf, "fasta"))
    else:
        with open(fh, mode = 'rt') as f:
            content = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    return(content)

def write_cds_faa(strain):
    """
    Takes the assembly.info file as a data.frame, and goes into each assembly directory (GCF something-something) to create a cds.faa file.
    The strain.cds.faa file is made by comparing the cds_from_genomic.fna and the protein.faa file.
    cds_from_genomic.fna = genes that form coding sequences, ATGC alphabet (fna)
    protein.faa = accessioned protein products annotated on the genome assembly (amino acid alphabet)
    strain.cds.faa = headers from cds_from_genomic, with protein sequences from protein.faa
    This file is not the same as translated.faa, as this file is merely a naive translation of cds_from_genomic.
    This function, from the genomes/ directory, goes into each assembly dir and reads the contents of the .gz files of interest, does what it needs, then goes back in /genomes.

    :param info_df: pandas dataframe from assembly.info file
    :param logfile: filehandle for the logfile, which will contain the numbers of CDSs, proteins, and CDSs that are proteins.
    """
    proteins = fasta_gz_to_dict(pattern = f"*protein.faa*")
    cds = fasta_gz_to_dict(pattern = f"*cds_from_genomic.fna*")
    cds_faa = {}
    lt2number = {}
    number2wp = {}
    for key, value in cds.items():
        gene_number = re.search(r"_(\d+)$", key)
        strain_gene = f"{strain}_{gene_number.group(1)}"
        ## this works for fungi at least, should check if it still works for bacteria from RefSeq with the WP_ notation
        prot_match = re.search(
            r"\[protein_id=([A-Z]*_?[0-9]*\.[0-9]{1})\]", value.description
        )
        lt_match = re.search(
            r"\[locus_tag=([A-Z,0-9]*_[A-Z,0-9]*)\]", value.description, re.I
        )
        # match = re.search(
        #     "\[protein_id=(WP_[0-9]*\.[0-9]{1})\]", value.description)
        
        if prot_match is not None:
            pi = prot_match.group(1)
            cds_faa[strain_gene] = SeqRecord(
                #seq=proteins[pi].seq, id=f"{pi} {value.id}", name=value.name, description=value.description)
                seq=proteins[pi].seq, id=f"{strain_gene}|{pi}", name=value.name, description=value.description)
            number2wp[strain_gene] = pi
        else:
            continue
        if lt_match is not None:
            lt2number[lt_match.group(1)] = strain_gene
    pi_from_cds_faa = re.findall(
        r"cds_([A-Z]*_?[0-9]*\.[0-9]{1})", str(cds.keys())
    )
    # wp_from_cds_faa = re.findall(
    #     "(WP_[0-9]*\.[0-9]{1})", str(cds_faa.keys()))
    n_cds = len(cds.keys())
    n_prot = len(proteins.keys()) # same as unique items in wp_from_cds_faa, meaning the additional sequences in latter are paralog proteins
    n_cds_prot = len(pi_from_cds_faa)
    with open(f"{strain}.cds.faa", 'w') as fasta:
        SeqIO.write(cds_faa.values(), fasta, "fasta")
    return n_cds, n_prot, n_cds_prot, lt2number, number2wp

def write_cds_gff(strain, lt2number, number2wp):
    """
    This function modifies the cds.gff file to change the ID to the gene_number, based on the locus tag
    """
    ## If gff is gzipped, unzip it
    gff = glob.glob(f"*genomic.gff*")[0]
    if gff.endswith(".gz"):
        os.system('gunzip ' + gff)
    gff = gff.replace(".gz", "")

    ## Parse through the gff file manually, and: 
    ### if the line is a gene, identify the locus-tag (in Name=), define {number} (from lt2number), then replace ID=gene-(regex) with ID=gene-{number}
    ### if the line is a CDS, identify the locus-tag (in Parent=gene-(regex)) and WP (from ID=cds-(regex)), define {number} (from lt2number), check correspondance with WP and number2wp(number), then replace ID=cds-(regex) with ID=cds-{number} and Parent=gene-(regex) with Parent=gene-{number}
    with open(gff, 'r') as gff_fh_r:
        lines = gff_fh_r.readlines()
    with open(f"{strain}.cds.gff", 'w') as gff_fh_w:
        for line in lines:
            if line.startswith("#"):
                gff_fh_w.write(line)
            else:
                parts = line.split("\t")
                match parts[2]:
                    case "gene":
                        lt = re.search(r"locus_tag=([A-Z,0-9]*_[A-Z,0-9]*)", parts[8], re.I)
                        if lt is not None:
                            number = lt2number.get(lt.group(1), None)
                            if number is not None:
                                parts[8] = re.sub(r"ID=(gene-[A-Z,0-9]*_[A-Z,0-9]*)", f"ID=gene-{number}", parts[8])
                                gff_fh_w.write("\t".join(parts) + "\n")
                    case "CDS":
                        lt = re.search(r"Parent=gene-([A-Z,0-9]*_[A-Z,0-9]*)", parts[8], re.I)
                        wp = re.search(r"ID=cds-([A-Z]*_?[0-9]*\.[0-9]{1})", parts[8], re.I)
                        if lt is not None and wp is not None:
                            number = lt2number.get(lt.group(1), None)
                            if number is not None:
                                wp_number = number2wp.get(number, None)
                                if wp_number is not None and wp_number == wp.group(1):
                                    parts[8] = re.sub(r"ID=cds-[A-Z]*_?[0-9]*\.[0-9]{1}", f"ID=cds-{number}", parts[8])
                                    parts[8] = re.sub(r"Parent=gene-[A-Z,0-9]*_[A-Z,0-9]*", f"Parent=gene-{number}", parts[8])
                                    parts[8] = re.sub(f"protein_id=({wp_number})", f"protein_id={number}|{wp_number}", parts[8])
                                    gff_fh_w.write("\t".join(parts) + "\n")
                    case _:
                        continue

def write_correspondance_table(lt2number, number2wp, strain):
    """
    This function writes a correspondance table between locus tags and gene numbers, as well as gene numbers and WP numbers.
    """
    ## Transform dictionaries to dataframes
    lt2number_df = pd.DataFrame(lt2number.items(), columns=['locus_tag', 'gene_number'])
    number2wp_df = pd.DataFrame(number2wp.items(), columns=['gene_number', 'protein_id'])

    ## Merge dataframes on gene_number
    correspondance_df = pd.merge(lt2number_df, number2wp_df, on='gene_number', how='outer')

    ## Write the table to a tsv file
    correspondance_df.to_csv(f"{strain}.correspondance_table.tsv", sep='\t', index=False)

if __name__ == '__main__':
    main()