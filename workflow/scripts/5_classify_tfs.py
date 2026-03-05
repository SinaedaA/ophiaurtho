import pandas as pd
from Bio import SeqIO
import re
import argparse
import os
import glob

def main():
    dbcan_dir, cog_dir, outdir, tf_families, perc_threshold = parsing_arguments()
    ## If outdir does not exist, create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ## Glob all IPR tsv files (results of InterProScan)
    tf_hmm_files = glob.glob(f"{dbcan_dir}/*/TF_hmm_results.tsv")

    ## Read each table and store in list
    dataframes = [read_hmm_tables(tf_hmm) for tf_hmm in tf_hmm_files]
    tf_tab = pd.concat(dataframes)
    ## Write this combined table to a file
    tf_tab.to_csv(f"{outdir}/TF_hmm_results_combined.tsv", sep = "\t", index = False)

    ## Compute statistics for each strain: how many of each TF family in each strain
    tf_stats = tf_tab.groupby(["Strain", "HMM_name"])["GeneID"].nunique().reset_index()
    tf_stats.columns = ["strain", "Family_name", "n_tfs"]

    ## Get cog2seqids data 
    cog_seqids = get_cog_seqids(cog_dir, outdir)
    
    ## Merge tf_tab and cog_seqids, based on GeneID
    tf_tab = tf_tab.merge(cog_seqids, on="GeneID", how="left")
    tf_tab.drop(columns=["HMM_length", "Target_Length", "i-evalue", "HMM_from", "HMM_to", "Target_from", "Target_to"], inplace=True)

    ## Replace "NaN" in the HMM_name column with "None"
    tf_tab["HMM_name"] = tf_tab["HMM_name"].fillna("None")

    ## Create a new table, with COGs and the associated HMM names
    cog2tf = tf_tab.groupby(["cog", "HMM_name"])["GeneID"].nunique().reset_index()
    cog2tf.columns = ["cog", "HMM_name", "n"]
    cog2tf["total"] = cog2tf.groupby("cog")["n"].transform("sum")
    cog2tf["percentage"] = (cog2tf["n"] / cog2tf["total"]) * 100
    # none are "None", which is interesting... 

    ## example: 9_116|WP_189729662.1 is annotated as "0043669", which comes from GO, but when I run Interproscan on it, it's identified as a MarR TF... And MarR is actually in the dbCAN results, so it's not like it's absent... 
    for family in tf_families:
        ## Look for how many proteins were identified in this family
        print(f"Processing family: {family}")
        family_count = tf_tab[tf_tab["HMM_name"] == family].shape[0]
        print(f"Number of proteins in family {family}: {family_count}")
        ## Make the output directory for this family
        family_dir = os.path.join(outdir, family)
        print(f"Creating directories: {family_dir}/cogs/")
        if not os.path.exists(f"{family_dir}/cogs/"):
            os.makedirs(f"{family_dir}/cogs/")
        ## Get the COGs that have > 70% of genes from this family
        cog_family = cog2tf[(cog2tf["HMM_name"] == family) & (cog2tf["percentage"] > 70)]
        ## If there are more than 0 COGs, link the file to the family directory
        if cog_family.empty:
            print(f"No COGs found for family {family}")
            continue
        else:
            print(f"Found {cog_family.shape[0]} COGs for family {family}")
            ## Write the COGs to a file in this directory
            cog_family.to_csv(os.path.join(family_dir, f"{family}_cogs.tsv"), sep="\t", index=False)
            for cog in cog_family["cog"].unique():
                ## Link the COG.fasta file to the family directory
                if not os.path.exists(f"{family_dir}/cogs/{cog}.fasta"):
                    os.system(f"ln -s {cog_dir}/{cog}.fasta {family_dir}/cogs/{cog}.fasta")

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='5_classify_tfs.py',
        description='Identifies TFs from the dbCAN predictions.')

    parser.add_argument('--dbcan_dir',
                        help="Path to the directory containing the dbCAN predictions for each strain.")
    parser.add_argument('--cog_dir', 
                        help="Path to the directory containing the Orthogroup files, from the proteinortho run.")
    parser.add_argument('--outdir', default = "./", 
                        help="Path to output directory.")
    parser.add_argument('--tf_families', required = True, type=str,
                        help="Comma-separated list of TF families to classify. Example: LacI,GntR,IclR,TetR. If not provided, only LacI and GntR will be used.")
    parser.add_argument('--perc_threshold', default = 70, type=int,
                        help="Minimum percentage of proteins in a COG that have to be identified as the specific family. Default is 70%.")

    args = vars(parser.parse_args())
    dbcan_dir = os.path.abspath(args["dbcan_dir"])
    cog_dir = os.path.abspath(args["cog_dir"])
    outdir = os.path.abspath(args["outdir"])
    tf_families = args["tf_families"].split(",")
    perc_threshold = args["perc_threshold"]

    return dbcan_dir, cog_dir, outdir, tf_families, perc_threshold

def read_hmm_tables(tf_hmm):
    """
    This function reads a dbCAN HMM prediction table and returns a pandas DataFrame.
    The table is expected to have the following columns:
    GeneID, ProtID, HMM_name, HMM_length, Target_Name, Target_Length, i-evalue, HMM_from, HMM_to, Target_from, Target_to, Coverage
    """
    try:
        hmm_tab = pd.read_csv(tf_hmm, sep="\t", header=0)
        hmm_tab.drop(columns=["HMM File Name"], inplace=True)
        # change the column names to match the expected format
        hmm_tab.columns = ["HMM_name", "HMM_length", "Target_Name", "Target_Length", "i-evalue", "HMM_from", "HMM_to", "Target_from", "Target_to", "Coverage"]
        # add strain column, based on Target_Name (2_105|WP_114242398.1 - anything before the first "_")
        hmm_tab["Strain"] = hmm_tab["Target_Name"].apply(lambda x: x.split("_")[0])
        # add GeneID and ProtID columns, based on Target_Name (GeneID=everything before "|", ProtID=everything after "|")
        hmm_tab["GeneID"] = hmm_tab["Target_Name"].apply(lambda x: x.split("|")[0])
        hmm_tab["ProtID"] = hmm_tab["Target_Name"].apply(lambda x: x.split("|")[1] if "|" in x else x)
        # return the final table
        return hmm_tab
    except Exception as e:
        print(f"Error reading {tf_hmm}: {e}")
        return pd.DataFrame()

def get_cog_seqids(cog_dir, outdir):
    """
    This function reads each Orthogroup.fasta file, and returns a list of seqids (gene_number) contained in each Orthogroup.
    """
    with open(f"{outdir}/cog2seqid.csv", mode = 'w') as out:
        out.write("cog;id\n")
    fastaList = glob.glob(f"{cog_dir}/*.fasta")
    for fasta in fastaList:
        cog = fasta.replace(".fasta", "").replace(f"{cog_dir}/", "")
        with open(fasta, mode = 'rt') as f:
            fastdic = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        for key in fastdic.keys():
            with open(f"{outdir}/cog2seqid.csv", mode = 'a') as out:
                out.write(f"{cog};{key}\n")
    # read the file and return a pandas DataFrame
    cog_seqids = pd.read_csv(f"{outdir}/cog2seqid.csv", sep = ";", header = 0)
    cog_seqids['GeneID'] = cog_seqids["id"].apply(lambda x: x.split("|")[0])
    cog_seqids['ProtID'] = cog_seqids["id"].apply(lambda x: x.split("|")[1] if "|" in x else x)
    # write new file
    cog_seqids.to_csv(f"{outdir}/cog2seqid.tsv", sep = "\t", index = False)
    return cog_seqids

def classify_tfs(nested_ipr_dict, tf_tab):
    """
    This function takes a the nested IPR dictionary (strain => geneid => [ipr_hits]) and a table of TF domains.
    It will loop over each element in the ipr_dict, and merge the tables based on the IPR domains.
    It returns a dictionary with geneids as keys and a list of TF domain names as values.
    """
    full_tf_tab = pd.DataFrame(columns = ["geneid", "ipr_id", "Family_name", "Domain_name", "strain"])
    ## Loop over each strain
    for strain, pot_tf_tab in nested_ipr_dict.items():
        print(f"Classifying TFs for strain {strain}")
        tab_merged = pot_tf_tab.join(tf_tab.set_index("InterPro"), on = "ipr_id", how = "left")
        tab_merged[["strain"]] = strain
        tab_merged.drop(columns = ["Family_description", "Pfam"], inplace = True)
        full_tf_tab = pd.concat([full_tf_tab, tab_merged], ignore_index=True)
    
    return full_tf_tab

def flatten(a_list):
    return [x for xs in a_list for x in xs]

if __name__ == '__main__':
    main()