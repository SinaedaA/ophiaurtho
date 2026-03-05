import argparse
import re
import pandas as pd
import os
import glob

def main():
    tf_data, outpath = parsing_arguments()
    # Make outpath directory if it doesn't exist
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # Get all meme input files for parsing
    meme_files = glob.glob(f"{tf_data}/*/meme/*/meme.txt")

    # Loop over the meme files and parse them
    for meme_file in meme_files:
        # Define base of the output files
        basename = meme_file.split('/')[-2]
        out_base = re.sub("meme_", "", basename)
        # Read meme.txt file and store in dictionary
        with open(meme_file, 'r') as f:
            # loop over the index of the lines in the file
            line_dict = {}
            for i, line in enumerate(f.readlines()):
                line_dict[i] = line.strip()
        
        # Call the function to parse the meme file
        parse_meme_file(line_dict, outpath, out_base)

def parsing_arguments():
    parser = argparse.ArgumentParser(description="Extract information from MEME.txt files.")
    parser.add_argument("--tf_data", required=True, help="Directory containing TF data.")
    parser.add_argument("--outpath", required=True, help="Output folder to write extracted MEME information.")

    args = vars(parser.parse_args())
    tf_data = os.path.abspath(args["tf_data"])
    outpath = args["outpath"]

    return tf_data, outpath

def parse_meme_file(line_dict, outpath, out_base):
    # Loop over the lines and extract the relevant information
    check = False
    for key, line in line_dict.items():
        if line.startswith('MOTIF') and check == False:
            check = True
            consensus = line.split()[1]
            motif_n = f"motif{line.split()[2].split('-')[1]}"
            e_val = line.split()[-1]
            width = line.split()[5]
            n_sites = int(line.split()[8])
            # define outfile names
            loc_outfile = f"{outpath}/{out_base}_{motif_n}.tsv"
            mat_outfile = f"{outpath}/{out_base}_{motif_n}.fasta"
            continue
        if check == True and re.match(r'\s*.* sites sorted by position p-value', line):
            location_keys = range(key + 4, key + 4 + n_sites)
            location_elements = [line_dict[location_key].split() for location_key in location_keys]
            # keep only elements 0 to 3
            location_elements = [element[:4] for element in location_elements]
            # list to dataframe
            df = pd.DataFrame(location_elements, columns = ['Sequence name', 'Strand', 'Start', 'P-value'])
            df['motif'] = motif_n
            df['consensus'] = consensus
            df['e_value'] = e_val
            df['width'] = width
            with open(loc_outfile, 'a') as location_f:
                df.to_csv(location_f, sep="\t", index=False, header=location_f.tell()==0)
            continue
        if check == True and re.match(r'\s*.* in BLOCKS format', line):
            matrix_keys = range(key + 3, key + 3 + int(n_sites))
            with open(mat_outfile, 'w') as matrix_f:
                for matrix_key in matrix_keys:
                    matrix_f.write(f">{line_dict[matrix_key].split()[0]}\n")
                    matrix_f.write(f"{line_dict[matrix_key].split()[-2]}\n")
            check = False
            continue

if __name__ == "__main__":
    main()