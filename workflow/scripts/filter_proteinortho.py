import pandas as pd
import re
import argparse
import os
import glob
from Bio import SeqIO

def main():
    main_tsv, filter_type = parsing_arguments()

    ## Read proteinortho main TSV file
    tsv = pd.read_csv(main_tsv, sep="\t", header=0, dtype=str)
    # add first column, as a rowId
    tsv.insert(0, 'rowId', [f"OG_{i}" for i in range(1, 1 + len(tsv))])
    # get all "cds.faa" column names
    cds_cols = [col for col in tsv.columns if re.search(r'cds\.faa', col)]
    print(cds_cols)
    # pivot longer the columns that contain the "cds.faa" pattern
    tsv = tsv.melt(id_vars=['rowId', '# Species', 'Genes', 'Alg.-Conn.'], value_vars = cds_cols, var_name='genome', value_name='protein')
    # then, divide the protein column into rows, for each comma-separated value inside each cell
    tsv = tsv.assign(protein=tsv['protein'].str.split(',')).explode('protein').reset_index(drop=True)
    # filter out rows where protein is NaN
    print(tsv[tsv['rowId'] == 'OG_2'])


def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='filter_proteinortho.py',
        description='Filters proteinortho results based on user-defined criteria.')

    parser.add_argument('--main-tsv',
                        help="Path to the main proteinortho results TSV file.")
    parser.add_argument('--filter',
                        help="What kind of filter to apply ('none', 'targets', 'family').")
    args = vars(parser.parse_args())
    main_tsv = os.path.abspath(args["main_tsv"])
    filter_type = args["filter"]

    return main_tsv, filter_type

if __name__ == '__main__':
    main()