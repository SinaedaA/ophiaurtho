import pandas as pd
import glob
from Bio import SeqIO
import re
import argparse
import os
import gffutils
import warnings

def main():
    dbcan_dir, gff_dir, resources_dir, tf_class, method, extension, neighbours, outdir = parsing_arguments()
    if method == "neighbouring" and extension != 0:
        extension = 0
        warnings.warn("Selected method is 'neighbouring', CGC limit extension set to 0.")
    
    # Create output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Read GFF database files
    reduced_cds_dict = read_gff_dict(gff_dir)

    # Read and merge: TF classification
    tf_data = pd.read_csv(tf_class, sep="\t")
    # Get TF gene coordinates from reduced_cds_dict
    tf_data['start'] = tf_data['geneid'].apply(lambda x: reduced_cds_dict[x]['start'])
    tf_data['end'] = tf_data['geneid'].apply(lambda x: reduced_cds_dict[x]['end'])
    tf_data['contig'] = tf_data['geneid'].apply(lambda x: reduced_cds_dict[x]['contig'])
    tf_data['strand'] = tf_data['geneid'].apply(lambda x: reduced_cds_dict[x]['strand'])
    tf_data['description'] = tf_data['geneid'].apply(lambda x: reduced_cds_dict[x]['description'])

    # Glob dbCAN output files
    cgc_files = glob.glob(f"{dbcan_dir}/*/cgc_standard_out.tsv")
    overview_files = glob.glob(f"{dbcan_dir}/*/overview.tsv")
    substrate_pred_files = glob.glob(f"{dbcan_dir}/*/substrate_prediction.tsv")
    
    # Concatenate overview files
    overview = concat_overview_file(overview_files)
    substr_pred = concat_substrate_pred_files(substrate_pred_files)

    # Read the fam-substrate-mapping.tsv
    fam_substrate_mapping = pd.read_csv(f"{resources_dir}/fam-substrate-mapping.tsv", sep="\t")
    fam_substrate_mapping.columns = fam_substrate_mapping.columns.str.strip()
    fam_substrate_mapping = fam_substrate_mapping[['Substrate_curated', 'Family', 'Name', 'EC_Number']]

    # Concatenate CGC files and create cgc_extended
    dataframes = [read_cgc_out(cgc_tsv, extension) for cgc_tsv in cgc_files] # extension is 0 when the method is neighbouring
    cgc_extended = pd.concat(dataframes, ignore_index=True)
    # Remove TFs from cgc_extended (as I don't trust the dbCAN assigned TFs)
    cgc_extended_no_tfs = cgc_extended[cgc_extended['Gene Type'] != 'TF']
    # Add TFs to the CGC clusters, based on overlap
    cgc_with_tfs = check_overlap(cgc_extended_no_tfs, tf_data)

    # Add TFs outside of the CGC limits, if method is neighbouring
    if method == "neighbouring" or method == "both":
        cgc_w_neighTFs = get_neighbouring(cgc_with_tfs, reduced_cds_dict, neighbours, tf_data)
        cgc_w_neighTFs.drop(columns=['CGC#', 'strain'], inplace=True)
        print(cgc_w_neighTFs[cgc_w_neighTFs['Gene Type'].isin(['Neighbouring TF', 'Overlapping TF'])]['Gene Type'].value_counts())
        cgc_new = cgc_w_neighTFs.copy()
    else:
        cgc_new = cgc_with_tfs.copy()

    # Merge fam-substrate-mapping with overview
    overview_substr = overview.merge(fam_substrate_mapping, left_on="Simplified_domain", right_on="Family", how="left")
    overview_substr = overview_substr[['ProteinID', 'Predicted_domain', 'Family', 'Substrate_curated', 'Name', 'EC_Number']]

    # Merge cgc_with_tfs with substr_pred
    cgc_substr = cgc_new.merge(substr_pred[['strain_CGC', 'Substrate', 'Substrate_score']], on=['strain_CGC'], how='left')
    
    # Merge cgc_final with overview_substr
    cgc_final = cgc_substr.merge(overview_substr[['ProteinID', 'Predicted_domain', 'Substrate_curated', 'Name', 'EC_Number']], left_on=['Protein ID'], right_on=['ProteinID'], how='left')
    cgc_final.drop(columns = ['ProteinID'], inplace=True)

    # Rename columns
    cgc_final.rename(columns={'Substrate': 'CGC_Substrate', 
                              'Substrate_score': 'CGC_Substrate_score',
                              'Substrate_curated': 'Gene_Substrate',
                              'Name': 'Name',
                              'EC_Number': 'EC_Number',
                              'Protein ID': 'ID'}, inplace=True)
    cgc_final = cgc_final[['ID', 'Contig ID', 'Gene Type', 'Gene Start', 'Gene Stop', 'Gene Strand', 'Gene Annotation', 'Predicted_domain', 'Name', 'Gene_Substrate', 'EC_Number',
                          'strain_CGC', 'CGC start', 'CGC stop', 'CGC_Substrate', 'CGC_Substrate_score', 'TF family']]
    
    # Save the final DataFrame to a tsv file
    output_file = os.path.join(outdir, "dbCAN_TFs_summary.tsv")
    cgc_final.to_csv(output_file, sep="\t", index=False)
    print(f"Output saved to {output_file}")

def parsing_arguments():
    """
    Parses the arguments from the command line.
    """
    parser = argparse.ArgumentParser(description="Integrate dbCAN and TF data.")
    parser.add_argument("--dbcan_dir", type=str, required=True, 
                        help="Path to the dbCAN output folder.")
    
    parser.add_argument("--db_dir", type=str, required=True, 
                        help="Path to the dbCAN db directory, containing fam-substrate-mapping.tsv.")
    
    parser.add_argument("--gff_db_dir", type = str, required=True,
                        help="Directory containing the strain.gff.db files, created during the make_cds_faa_gff step. Should be same directory as the strain.faa and strain.gff files.")
    
    parser.add_argument("--tf_classification", type=str, required=True, 
                        help="Path to the TF classification file.")
    
    parser.add_argument("--method", type=str, default="overlap", choices=["overlap", "neighbouring", "both"], 
                        help="Method to check if a TF should be added into a CGC cluster. Default is 'overlap'.")
    
    parser.add_argument("--extension", type = int, default = 100, 
                        help = "Number of bp to extend the CGC cluster limits by, if method = overlap. " \
                        "Default is 100 bp. Ignored if method is neighbouring.")
    
    parser.add_argument("--neighbours", type = int, default = 1, 
                        help = "Number of neighbouring genes to consider if method = neighbouring. " \
                        "At default=1, we only look at the previous and the next gene (compared to the CGC members). " \
                        "Ignored if method is overlap.")
    
    parser.add_argument("--outdir", required = True, 
                        help = "Path to the output directory")
    
    args = parser.parse_args()
    dbcan_dir = os.path.abspath(args.dbcan_dir)
    resources_dir = os.path.abspath(args.db_dir)
    gff_dir = os.path.abspath(args.gff_db_dir)
    tf_class = os.path.abspath(args.tf_classification)
    method = args.method
    extension = args.extension
    neighbours = args.neighbours
    outdir = os.path.abspath(args.outdir)

    return dbcan_dir, gff_dir, resources_dir, tf_class, method, extension, neighbours, outdir

def read_gff_dict(gff_dir):
    """
    This function takes a directory containing the strain.gff.db files (created during "make_cds_faa_gff" step). 
    1. Create gff_dict and gff_cds_dict, to further create reduced dictionaries. 
        - gff_dict: Keys = strains, values = db (features = gene) (gffutils)
        - gff_cds_dict: Keys = strains, values = db (features = CDS) (gffutils)
    2. Reduced dictionaries for genes and CDSs are created and returned
        - reduced_gene_dict : keys = locus_tag, values = dict(contig, start, end, strand))
        - reduced_cds_dict : keys = accession, values = dict(contig, start, end)
    """
    gff_dict = {}
    gff_cds_dict = {}

    print(f"Looking for GFF databases in {gff_dir}")
    gff_dbs = glob.glob(f"{gff_dir}/*.gff.db")
    for gff_db in gff_dbs:
        key = gff_db.split("/")[-1].split(".")[0]
        gff = gffutils.FeatureDB(gff_db, keep_order = True)
        gff_dict[key] = gff.all_features(featuretype="gene")
        gff_cds_dict[key] = gff.all_features(featuretype="CDS")
            #print(gff.all_features(featuretype = "CDS")) ## generator
    ## create a reduced_cds_dict as well, with accession -> contig, start, end
    reduced_cds_dict = {}
    for generator in gff_cds_dict.values():
        for feat in generator:
            contig = feat.seqid
            #accession = feat.id.split("-")[1] # remove the cds- part
            protid = feat.attributes['protein_id'][0]
            start = feat.start
            end = feat.end
            strand = feat.strand
            description = feat.attributes.get('product', [''])[0]  # Get product description if available
            reduced_cds_dict[protid] = {"contig": contig, "start": start, "end": end, "strand": strand, "description": description}
    return reduced_cds_dict

def read_cgc_out(cgc_file, extension):
    """
    Reads a tsv file using pandas, and returns a pd.DataFrame object. 
    For each CGC_standard_out file, it reads it and add a column for CGC start and stop coordinates, + or - the extension parameter (integer).
    I also add a column for the strain_CGC ID, as otherwise we have the same CGC# in different strains, and we end up concatenating all the tables. 
    Returns the modified CGC DataFrame (pandas).
    """
    cgc = pd.read_csv(cgc_file, sep = "\t")
    cgc['strain'] = cgc['Protein ID'].apply(lambda x: x.split("_")[0])
    cgc['strain_CGC'] = cgc["strain"] +"_"+ cgc["CGC#"].astype(str)
    cgc['CGC start'] = cgc.groupby(['strain_CGC'])['Gene Start'].transform('min') - extension
    cgc['CGC stop'] = cgc.groupby(['strain_CGC'])['Gene Stop'].transform('max') + extension
    return(cgc)

def concat_overview_file(overview_files):
    """
    Reads an overview file and returns a DataFrame with the relevant columns.
    The overview file contains information about CGC clusters, including their start and stop coordinates.
    """
    dataframes = [pd.read_csv(overview, sep = "\t") for overview in overview_files]
    overview = pd.concat(dataframes, ignore_index=True)
    # Add a strain column
    overview['strain'] = overview['Gene ID'].apply(lambda x: x.split("_")[0])
    # Drop all columns, except strain, Gene ID, Recommended Results
    overview = overview[['strain', 'Gene ID', 'Recommend Results']]
    # Rename columns
    overview.columns = ['strain', 'ProteinID', 'Predicted_domain']
    # Flatten the 'Predicted_domain' column, to have one row per domain
    overview['Predicted_domain'] = overview['Predicted_domain'].str.split('|')
    overview = overview.explode('Predicted_domain')
    # Make a simplified version of the 'Predicted_domain' column
    overview['Simplified_domain'] = overview['Predicted_domain'].apply(lambda x: x.split("_")[0] if isinstance(x, str) else x)
    return overview

def concat_substrate_pred_files(substrate_pred_files):
    """
    Reads substrate prediction files and returns a DataFrame with the relevant columns.
    The substrate prediction files contain information about the predicted substrates for each CGC cluster.
    """
    dataframes = [add_strain_column(prediction) for prediction in substrate_pred_files]
    substr_pred = pd.concat(dataframes, ignore_index=True)

    substr_pred['CGC'] = substr_pred['#cgcid'].apply(lambda x: x.split("|")[1])
    substr_pred['strain_CGC'] = substr_pred['strain'] + "_" + substr_pred['CGC']
    
    substr_pred = substr_pred[['strain', 'strain_CGC', 'dbCAN-sub substrate', 'dbCAN-sub substrate score']]
    substr_pred.columns = ['strain', 'strain_CGC', 'Substrate', 'Substrate_score']

    return substr_pred

def add_strain_column(path):
    strain = os.path.basename(os.path.dirname(path))
    df = pd.read_csv(path, sep="\t")
    df['strain'] = strain
    return df

def check_overlap(cgc_df, tf_df):
    """
    Checks for overlap between CGC cluster limits and TF coordinates.
    Returns a pd.DataFrame with overlaps.
    """
    ## Create a simplified version, containing only the strain, strain_CGC, has_tf, CGC start and CGC stop columns
    cgc_simple = cgc_df[['strain', 'strain_CGC', 'Contig ID', 'CGC start', 'CGC stop']].drop_duplicates()
    
    ## Print some data: number of unique CGC clusters, number of CGC clusters with dbCAN detected TFs
    print(f"Found {len(cgc_simple)} CGC clusters across all strains.") 

    ## Order both tables by coordinates, per strain
    cgc_simple = cgc_simple.sort_values(by=['strain', 'CGC start'])
    tf_df = tf_df.sort_values(by=['strain', 'start'])

    ## Create dictionaries for each table
    cgc_dict = make_cgc_dict(cgc_simple)
    tf_dict = make_tf_dict(tf_df)

    ## Check for overlaps between extended CGC clusters and TF coordinates
    for strain, tfs in tf_dict.items():
        if strain not in cgc_dict:
            continue
        cgc_strain = cgc_dict[strain]
        for tf in tfs:
            for cgc_id, cgc_info in cgc_strain.items():
                ## could replace this with ranges ???
                if (tf['start'] <= cgc_info['stop'] and tf['stop'] >= cgc_info['start']):
                    # add the TF to the CGC DF as well
                    cgc_df = pd.concat([cgc_df, pd.DataFrame({
                        'CGC#': cgc_id.split("_")[1],  # Extract CGC number from strain_CGC
                        'Gene Type': 'Overlapping TF',
                        'Contig ID': cgc_info['contig'],  # Assuming the contig ID is part of the geneid
                        'Protein ID': [tf['geneid']],
                        'Gene Start': [tf['start']],
                        'Gene Stop': [tf['stop']],
                        'Gene Strand': [tf['strand']],  # Assuming strand is always '+', can be modified if needed
                        'Gene Annotation': [tf['description']],  # Placeholder, can be filled with actual annotation if available
                        'strain': strain,
                        'strain_CGC': cgc_id,
                        'CGC start': [cgc_info['start']],
                        'CGC stop': [cgc_info['stop']],
                        'TF family': [tf['family']]  # Add TF family if available
                    })])
    ## Sort the cgc_df by strain_CGC and Gene Start
    cgc_df = cgc_df.sort_values(by=['strain_CGC', 'Gene Start']).reset_index(drop=True)
    return cgc_df

def get_neighbouring(cgc_df, reduced_cds_dict, neighbours, tf_data):
    '''
    Checks for neighbouring genes of CGC clusters.
    Returns the same CGC DataFrame, but with added neighbouring genes for each CGC cluster, and updated CGC coordinates.
    '''
    # cgc_df has "Protein ID" column, that corresponds to the keys reduced_cds_dict
    ki = dict()
    ik = dict()
    for i, k in enumerate(reduced_cds_dict):
        ki[k] = i # k = keys of reduced_cds_dict, values is the index
        ik[i] = k # i = index, values = keys of reduced_cds_dict
    # 1. For each CGC, get the Protein ID with the smallest coordinates (FIRST) and highest coordinates (LAST), i.o.w. the genes delimiting the cluster
    ## filter cgc_df based on start and stop (must be equal to CGC start and stop)
    cgc_firsts_lasts = cgc_df[(cgc_df['Gene Start'] == cgc_df['CGC start']) | (cgc_df['Gene Stop'] == cgc_df['CGC stop'])][['strain_CGC', 'Protein ID', 'Gene Start', 'Gene Stop', 'CGC start', 'CGC stop']]
    cgc_firsts_lasts['Position'] = cgc_firsts_lasts.apply(lambda row: 'FIRST' if row['Gene Start'] == row['CGC start'] else 'LAST', axis=1)
    ## Get all unique strain_CGC IDs
    unique_strain_cgc = cgc_firsts_lasts['strain_CGC'].unique()
    ## Loop over them, and for each strain_CGC, get the FIRST and LAST Protein IDs, and find the one before FIRST and the one after LAST
    for strain_cgc in unique_strain_cgc:
        strain = strain_cgc.split("_")[0] # get the strain name from the strain_CGC ID
        first_last_df = cgc_firsts_lasts[cgc_firsts_lasts['strain_CGC'] == strain_cgc]
        if len(first_last_df) == 2: # should always be true
            first = first_last_df[first_last_df['Position'] == 'FIRST']['Protein ID'].values[0]
            last = first_last_df[first_last_df['Position'] == 'LAST']['Protein ID'].values[0]
            first_i = ki[first]
            last_i = ki[last]

            # 2. Get the Protein IDs of the neighbours, i.e. the X ones before FIRST and after LAST
            before_first = [ik[first_i - n] for n in range(1, neighbours + 1) if first_i - n >= 0]
            after_last = [ik[last_i + n] for n in range(1, neighbours + 1) if last_i + n < len(ik)] 

            # 3. Check if any of the neighbours are TFs
            tf_neighbours = tf_data[tf_data['geneid'].isin(before_first + after_last)]['geneid'].tolist()
            
            if len(tf_neighbours) > 0:
                # 3. Create a DataFrame with the neighbouring TFs only
                cgc_extra_genes = pd.DataFrame({
                    'CGC#': strain_cgc.split("_")[1], # Extract CGC number from strain_CGC
                    'Gene Type': 'Neighbouring TF', # single value, repeated
                    'Contig ID': reduced_cds_dict[first]['contig'], # single value, repeated
                    'Protein ID': tf_neighbours, # list of Protein IDs
                    'Gene Start': [reduced_cds_dict[before_pid]['start'] for before_pid in tf_neighbours],
                    'Gene Stop': [reduced_cds_dict[before_pid]['end'] for before_pid in tf_neighbours],
                    'Gene Strand': [reduced_cds_dict[before_pid]['strand'] for before_pid in tf_neighbours],
                    'Gene Annotation': [reduced_cds_dict[before_pid]['description'] for before_pid in tf_neighbours],
                    'strain': strain,
                    'strain_CGC': strain_cgc,
                    'CGC start': cgc_firsts_lasts['CGC start'].values[0], # single value, repeated
                    'CGC stop': cgc_firsts_lasts['CGC stop'].values[0], 
                    'TF family': tf_data[tf_data['geneid'].isin(tf_neighbours)]['Short_name'].tolist()
                })
                # 4. Append the neighbours to the cgc_df
                cgc_df = pd.concat([cgc_df, cgc_extra_genes], ignore_index=True)
    # 5. Sort the cgc_df by strain_CGC and Gene Start
    cgc_df = cgc_df.sort_values(by=['strain_CGC', 'Gene Start']).reset_index(drop=True)
    # 6. Group the DF by strain_CGC, and change the CGC start and stop coordinates to the min and max of the Gene Start and Gene Stop coordinates
    cgc_df['CGC start'] = cgc_df.groupby('strain_CGC')['Gene Start'].transform('min')
    cgc_df['CGC stop'] = cgc_df.groupby('strain_CGC')['Gene Stop'].transform('max')

    # Return the updated cgc_df
    return cgc_df
        
def make_cgc_dict(cgc_df):
    """
    Creates a nested dictionary from the CGC DataFrame.
    First level keys are strains, and the second level keys are strain_CGC IDs.
    Each value is a dictionary with 'start' and 'stop' coordinates for the CGC cluster.
    Returns the nested dictionary.
    """
    cgc_dict = {}
    for _, row in cgc_df.iterrows():
        strain = row['strain']
        cgc_id = row['strain_CGC']
        if strain not in cgc_dict:
            cgc_dict[strain] = {}
        cgc_dict[strain][cgc_id] = {'start': row['CGC start'], 'stop': row['CGC stop'], 'contig': row['Contig ID']}
    return cgc_dict

def make_tf_dict(tf_df):
    """
    Creates a dictionary from the TF DataFrame.
    First level keys are strains, and each value is a list of dictionaries with 'geneid', 'start', 'stop', and 'family' for each TF.
    Returns the dictionary.
    """
    tf_dict = {}
    for _, row in tf_df.iterrows():
        strain = row['geneid'].split("_")[0]
        geneid = row['geneid']
        if strain not in tf_dict:
            tf_dict[strain] = []
        # Append to as a list of list
        tf_dict[strain].append({
            'geneid': geneid,
            'start': min(row['start'], row['end']),
            'stop': max(row['start'], row['end']),
            'strand': row['strand'],
            'family': row['Short_name'],
            'description': row['description']
        })
    return tf_dict

def range_intersect(r1, r2):
    return range(max(r1.start,r2.start), min(r1.stop,r2.stop)) or None


if __name__ == '__main__':
    main()