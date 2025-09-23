import pandas as pd
import re
import argparse
import os
import glob
from Bio import SeqIO

def main():
    ipr_dir, tf_domains, outdir, proteinortho_tsv, tf_families, perc_threshold = parsing_arguments()
    ## If outdir does not exist, create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    ##### INTERPRO DATA #####
    ## Glob all IPR tsv files (results of InterProScan)
    ipr_files = glob.glob(f"{ipr_dir}/*.tsv")

    ##### TF DOMAIN DATA #####
    ## Make TF table and TF IPR domain list
    tf_tab = pd.read_csv(tf_domains, sep = "\t", header = 0)
    tf_ipr_list = [item.split(',') for item in list(tf_tab.InterPro)]
    tf_ipr_list = flatten(tf_ipr_list)

    ###### CLASSIFYING TFs ######
    ## Make full dictionary for each strain, and each geneid, and all associated IPR hits (strain => geneid => IPR hits) (nested dict)
    ipr_full_dict = make_strain_ipr_dict(ipr_files, tf_ipr_list)
    tf_full_tab = classify_tfs(ipr_full_dict, tf_tab)
    
    ## Return stats for each strain and each family of TFs
    tf_stats = tf_full_tab.groupby(["strain", "Short_name"])["geneid"].nunique().reset_index()
    tf_stats.columns = ["strain", "Short_name", "n_tfs"]
    
    ## Write both tables to file
    tf_full_tab.to_csv(f"{outdir}/TF_classification.tsv", sep = "\t", index = False)
    tf_stats.to_csv(f"{outdir}/TF_stats.tsv", sep = "\t", index = False)   

    ###### COG DATA ######
    ## Get the COG to seqids mapping file, based on proteinortho.tsv
    proteinortho_tsv, cog2seqid = make_cog2seqid(proteinortho_tsv)

    ## Merge tf_tab and cog_seqids, based on GeneID
    tf_full_tab = tf_full_tab.merge(cog2seqid, on="geneid", how="left")
    
    ## Create a new table, with COGs and the associated HMM names
    cog2tf = tf_full_tab.groupby(["cog", "Short_name", "n_genes"])["geneid"].nunique().reset_index()
    cog2tf.columns = ["cog", "Short_name", "n_genes", "n"]
    # Make the n_genes column numeric
    cog2tf["n_genes"] = pd.to_numeric(cog2tf["n_genes"])
    cog2tf["percentage"] = (cog2tf["n"] / cog2tf["n_genes"]) * 100
    ## Add this info to the tf_full_tab dataframe
    tf_full_tab = tf_full_tab.merge(cog2tf[["cog", "Short_name", "n", "percentage"]], on=["cog", "Short_name"], how="left")

    for family in tf_families:
        ## Look for how many proteins were identified in this family
        print(f"Processing family: {family}")
        ## Identify TFs for further analysis (percentage > threshold), and from the family of interest
        tfs_oi = tf_full_tab[(tf_full_tab["Short_name"] == family) & (tf_full_tab["percentage"] > perc_threshold)]
        cogs_oi = tfs_oi["cog"].unique()

        ## Check content of cogs_oi
        if cogs_oi.size == 0:
            print(f"No COGs found for family {family}")
            continue
        else:
            print(f"Found {len(cogs_oi)} COGs for family {family}")
            ## Get corresponding rowIds from the original proteinortho tsv file
            rowIds_oi = [re.sub('OG_', '', cog) for cog in cogs_oi]
            rowIds_oi = [int(rowId) - 1 for rowId in rowIds_oi] # -1 because python is 0-indexed

            ## Filter proteinortho_tsv, to only keep COGs that are identified as TFs (keep only rowIds_oi as row numbers)
            tf_tsv = proteinortho_tsv.iloc[rowIds_oi]
            tf_tsv.to_csv(f"{outdir}/{family}.proteinortho.tsv", sep = "\t", index = False)

            ## I wanted to write the COG fasta files, but for that I need the proteinortho conda environment.
    ## Write a file to say that the script finished without error
    with open(f"{outdir}/../.classify_tfs_done", mode = 'w') as f:
        f.write("done\n")

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='5_classify_tfs_IPR.py',
        description='Classifies transcription factors (TFs) based on their InterProScan prediction tables.')

    parser.add_argument('--ipr_dir',
                        help="Path to the directory containing InterProScan predictions.")
    parser.add_argument('--tf_domains',
                        help="Path to the tsv file containing the TF domains (IPR and Pfam).")
    parser.add_argument('--proteinortho_tsv',
                        help = "Path to the [project_name].proteinortho.tsv file.")
    parser.add_argument('--tf_families', required = True, type=str,
                        help="Comma-separated list of TF families to classify. Example: LacI,GntR,IclR,TetR. If not provided, only LacI and GntR will be used.")
    parser.add_argument('--perc_threshold', default = 70, type=int,
                        help="Minimum percentage of proteins in a COG that have to be identified as the specific family. Default is 70%.")
    parser.add_argument('--outdir', default = "./", 
                        help="Path to output directory.")

    args = vars(parser.parse_args())
    ipr_dir = os.path.abspath(args["ipr_dir"])
    tf_domains = os.path.abspath(args["tf_domains"])
    proteinortho_tsv = os.path.abspath(args["proteinortho_tsv"])
    tf_families = args["tf_families"].split(",")
    perc_threshold = args["perc_threshold"]
    outdir = os.path.abspath(args["outdir"])

    return ipr_dir, tf_domains, outdir, proteinortho_tsv, tf_families, perc_threshold

def ipr_to_dict(ipr_file, tf_ipr_list):
    """
    This function takes a filehandle of an InterProScan prediction table, then reads the content and filters the table. 
    Filtering happens like this: 
        1. Remove any hits where evalue <= 1e-5
        2. It only keeps the top 3 IPR hits.
        3. Keep only IPR domains associated with TFs.
    Then it stores the geneid and the corresponding top 3 IPR hits in a dictionary, which it returns.
    """
    ## Get the strain name from the file name
    strain = os.path.basename(ipr_file).split(".")[0]

    ## Read the table and set the column names
    ipr_tab = pd.read_csv(ipr_file, sep = "\t", header = None)
    ipr_tab.columns = ["geneid", "hash", "prot_length", "tool", "id", "description", "start", "end", "evalue", "status", "date", "ipr_id", "ipr_description", "na1", "na2"]
    
    ## 1. Filter the table for IPR domains (not "-") and for status == True (meaning the match has been integrated in InterPro) ##evalue <= 1e-5
    sub_ipr_tab = ipr_tab[["geneid", "ipr_id", "status", "evalue"]]
    filt_ipr_tab = sub_ipr_tab[sub_ipr_tab.ipr_id != "-"]
    filt_ipr_tab = filt_ipr_tab[filt_ipr_tab.status == "T"]
    ipr_hit_n = len(filt_ipr_tab)
    print(f"Strain {strain} has {ipr_hit_n} IPR domain hits for {len(filt_ipr_tab.geneid.unique())} unique geneids (this is NOT the total number of genes in the genome).")
    ## 2. Keep only top 3 IPR domains
    filt_ipr_tab = filt_ipr_tab.groupby("geneid",as_index=False)[['geneid', 'ipr_id', 'evalue']].apply(lambda x: x.nsmallest(3, 'evalue'), include_groups=False)
    filt_ipr_tab.reset_index(inplace=True)
    filt_ipr_tab.drop(columns=["level_0","level_1"], inplace=True)
    print(f"After filtering, we kept {len(filt_ipr_tab)}/{ipr_hit_n} IPR domain hits for {len(filt_ipr_tab.geneid.unique())} unique geneids.")
    ## 3. Keep only proteins containing IPR domains associated with TFs
    potential_tfs = filt_ipr_tab[filt_ipr_tab.ipr_id.isin(tf_ipr_list)].geneid.unique()
    print(f"Of those, {len(potential_tfs)} contain a DNA-binding domain associated with a TF.")
    pot_tf_tab = filt_ipr_tab[filt_ipr_tab.ipr_id.isin(tf_ipr_list)][["geneid", "ipr_id"]].drop_duplicates() ## issue with genes that have 2 DBDs

    return pot_tf_tab

def make_strain_ipr_dict(ipr_files, tf_ipr_list):
    """
    This function takes a list of InterProScan prediction tables, then reads the content and stores it in a dictionary, which it returns.
    The keys are the strain names (GCF_XXXX), and the values are dictionaries with geneids as keys and a list of IPR values as values.
    """
    strain_tf_dict = {}
    for ipr_file in ipr_files:
        strain = os.path.basename(ipr_file).split(".")[0]
        pot_tf_tab = ipr_to_dict(ipr_file, tf_ipr_list)
        strain_tf_dict[strain] = pot_tf_tab

    return strain_tf_dict

def classify_tfs(nested_ipr_dict, tf_tab):
    """
    This function takes a the nested IPR dictionary (strain => geneid => [ipr_hits]) and a table of TF domains.
    It will loop over each element in the ipr_dict, and merge the tables based on the IPR domains.
    It returns a dictionary with geneids as keys and a list of TF domain names as values.
    """
    full_tf_tab = pd.DataFrame(columns = ["geneid", "ipr_id", "Family_name", "Short_name", "Domain_name", "strain"])
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

# def get_cog_seqids(cog_dir, outdir):
#     """
#     This function reads each Orthogroup.fasta file, and returns a list of seqids (gene_number) contained in each Orthogroup.
#     """
#     with open(f"{outdir}/cog2seqid.csv", mode = 'w') as out:
#         out.write("cog;id;start;end;strand\n")
#     fastaList = glob.glob(f"{cog_dir}/*.fasta")
#     for fasta in fastaList:
#         cog = fasta.replace(".fasta", "").replace(f"{cog_dir}/", "")
#         with open(fasta, mode = 'rt') as f:
#             fastdic = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
#         for key, value in fastdic.items():
#             description = value.description
#             location = re.search(r"location=(.*)\]{1}? ", description)[1]
#             if re.search(r"<|>", location): continue
#             start, end, strand = parse_location(location)
#             with open(f"{outdir}/cog2seqid.csv", mode = 'a') as out:
#                 out.write(f"{cog};{key};{start};{end};{strand}\n")
#     # read the file and return a pandas DataFrame
#     cog_seqids = pd.read_csv(f"{outdir}/cog2seqid.csv", sep = ";", header = 0)
#     cog_seqids['GeneNumber'] = cog_seqids["id"].apply(lambda x: x.split("|")[0])
#     cog_seqids['ProtID'] = cog_seqids["id"].apply(lambda x: x.split("|")[1] if "|" in x else x)
#     # Drop the start, end and strand columns, as they are not needed for the final output
#     cog_seqids.drop(columns=["start", "end", "strand"], inplace=True)
#     # write new file
#     cog_seqids.columns = ["cog", "geneid", "GeneID", "ProtID"]
#     cog_seqids.to_csv(f"{outdir}/cog2seqid.tsv", sep = "\t", index = False)
#     return cog_seqids

def make_cog2seqid(proteinortho_tsv):
    ## Read proteinortho main TSV file
    original_tsv = pd.read_csv(proteinortho_tsv, sep="\t", header=0, dtype=str)
    tsv = original_tsv.copy()
    # add first column, as a rowId
    tsv.insert(0, 'OG', [f"OG_{i}" for i in range(1, 1 + len(tsv))])
    # get all "cds.faa" column names
    cds_cols = [col for col in tsv.columns if re.search(r'cds\.faa', col)]
    # pivot longer the columns that contain the "cds.faa" pattern
    tsv = tsv.melt(id_vars=['OG', '# Species', 'Genes', 'Alg.-Conn.'], value_vars = cds_cols, var_name='proteome', value_name='protein')
    # then, divide the protein column into rows, for each comma-separated value inside each cell
    tsv = tsv.assign(protein=tsv['protein'].str.split(',')).explode('protein').reset_index(drop=True)
    # add a strain column, based on the proteome column (0.cds.faa --> 0)
    tsv['strain'] = tsv['proteome'].apply(lambda x: re.sub(r'\.cds\.faa$', '', x))
    # change colnames
    tsv.columns = ["cog", "n_species", "n_genes", "alg_conn", "proteome", "geneid", "strain"]

    return(original_tsv, tsv)

# def parse_location(location):
#     """
#     parse_location parses the output to the regex match on a cds.faa header, to extract the location information. 
#     This can look like: 
#     On positive strand: normal [11178..11507] || join(11854..12282,12340..12810)
#     On negative strand: [complement(876..1205)] || complement(join(3015..3209,3255..3981,4032..5107))
#     If the location contains "<" or ">", means it's a partial gene. Skipping those for now.
#     """
#     coords_str = re.findall(r'\d+', location) # returns string
#     coords = list(map(int, coords_str))
#     #print(coords) # [181763, 181846, 181897, 182256, 182311, 182691] -> group by 2 and then get range(2nd+1,3rd-1) etc
#     if re.search(r"^\d+\.\.\d+$|^join\(\d+.*\)", location):
#         strand = "+"
#         start = min(coords)
#         end = max(coords)
#     if re.search(r"^complement\(\d+\.\.\d+\)|^complement\(join\(\d+.+\)", location):
#         strand = "-"
#         start = max(coords)
#         end = min(coords)
#     return start, end, strand

if __name__ == '__main__':
    main()