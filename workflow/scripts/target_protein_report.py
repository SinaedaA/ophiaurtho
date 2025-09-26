import re
import argparse
import glob
import pandas as pd
import os
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import BytesIO
from fpdf import FPDF
from fpdf.fonts import FontFace
from fpdf.enums import TableCellFillMode

class PDF(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 16)
        self.set_fill_color(70, 80, 99) #dark grey/blue
        self.set_draw_color(0, 0, 0) #black
        self.set_text_color(250, 251, 252) #almost white
        self.cell(0, 10, self.title, new_x="LMARGIN", new_y="NEXT", align="L", fill = True)
        self.ln(3)
    def footer(self):
        # Position cursor at 1.5 cm from bottom:
        self.set_y(-15)
        # Setting font: helvetica italic 10
        self.set_font("helvetica", style="I", size=10)
        # Printing page number:
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="C")
    def subtitle1(self, text):
        self.set_font("Helvetica", size = 16, style = 'BU')
        self.set_text_color(0, 0, 0)
        self.cell(20, 10, text, new_x="LMARGIN", new_y="NEXT", align = 'L')
        self.ln(2)
    def subtitle2(self, text):
        self.set_font("Helvetica", size = 14, style = 'B')
        self.set_text_color(0, 0, 0) #black
        self.cell(20, 5, text, new_x="LMARGIN", new_y="NEXT", align = 'L')
        self.ln(2)
    def meme_section(self, text):
        self.set_font("Helvetica", "B", 14)
        self.set_fill_color(171, 194, 237) #dark grey/blue
        self.set_draw_color(0, 0, 0) #black
        self.set_text_color(250, 251, 252) #almost white
        self.cell(0, 7, text, new_x="LMARGIN", new_y="NEXT", align="C", fill = True)
        self.ln(1)
    def normal_text(self, text):
        self.set_font("Helvetica", size = 12)
        self.set_text_color(0, 0, 0) #black
        self.multi_cell(0, 5, text)
        self.ln()
    # def table(self, data, col_widths=None, row_height=10, align='C'):
    #     self.set_font("Helvetica", size = 12)
    #     self.set_text_color(0, 0, 0) #black
    #     for row in data:
    #         for i, item in enumerate(row):
    #             if col_widths:
    #                 self.cell(col_widths[i], row_height, str(item), border=1, align=align)
    #             else:
    #                 self.cell(40, row_height, str(item), border=1, align=align)
    #         self.ln(row_height)

def main():
    targets, tf_data, meme_info, outdir = parsing_arguments()
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Make a cog2seqid for the COGs that were written for the TF families of interest
    tf_cog2seqid = get_cog_seqids(tf_data)

    # Start PDF
    pdf = PDF()
    pdf.set_title("OphiAURTHO Report for Target Protein(s)")
    pdf.set_auto_page_break(auto = True, margin = 15)
    pdf.add_page()

    # Find the target proteins in the cog2seqid data
    target_cog2seqid = tf_cog2seqid[tf_cog2seqid["ProtID"].isin(targets)]
    target_cog2seqid_list = [target_cog2seqid.columns.values.tolist()] + target_cog2seqid.values.tolist()
    
    # Get the MEME results
    cogs = target_cog2seqid["cog"].unique().tolist()
    cog2png_df = get_cog2png_df(tf_data, meme_info, cogs)

    # Print the resulting table to the pdf report
    pdf.subtitle1("1. Target protein to Orthogroup mapping:")
    pdf.normal_text("The following section shows you which COG each of your target transcription factors are found in, along with the gene ID and protein ID.")
    # table
    pdf.set_draw_color(0, 0, 0)
    pdf.set_line_width(0.3)
    headings_style = FontFace(emphasis="BOLD", color=(250, 251, 252), fill_color=(70, 80, 99))
    with pdf.table(
        borders_layout="NO_HORIZONTAL_LINES",
        cell_fill_color=(224, 235, 255),
        cell_fill_mode=TableCellFillMode.ROWS,
        headings_style=headings_style,
        line_height=6,
        text_align=("CENTER", "CENTER", "CENTER", "CENTER"),
        width=180
    ) as table:
        for data_row in target_cog2seqid_list:
            row = table.row()
            for item in data_row:
                row.cell(str(item))
    pdf.write(5,"\n\n")
    
    i = 2
    # Make a subsection for each target protein
    for target in targets:
        ## Define COG and set new Header for each page
        cog = target_cog2seqid[target_cog2seqid["ProtID"] == target]["cog"].tolist()
        cog = cog[0] if len(cog) > 0 else None
        pdf.set_title(f"{target} | {cog}")
        pdf.add_page()
        if cog is None:
            pdf.normal_text(f"No OrthoGroup found for target protein {target}, skipping to next target.\nIt's possible that your target protein belongs to an OrthoGroup with less proteins than the min_prot parameter, set during the OphiAurtho pipeline. If this is the case, you can either adjust this parameter (not recommended to go below 4 proteins), or increase the number of genomes included in the analysis.")
            continue
        ###### SUBTITLE PER TARGET ######
        pdf.subtitle1(f"{i}. Results summary for {target} ({cog if cog is not None else 'None'})")
        ###### SUB-SUBTITLE -- PROTEIN DESCRIPTIONS ######
        pdf.subtitle2(f"{i}.1 Protein descriptions in the OrthoGroup ({cog})")
        #### Description ####
        pdf.normal_text("The following graph shows the distribution of protein descriptions for proteins, belonging to the same OrthoGroup as your target.")
        #### Plot Protein Descriptions ####
        cog_descriptions = get_protein_descriptions(cog, tf_data)
        plt.figure(figsize=(8, 6)) # create figure object
        sns.set_style("white")
        sns.set_style("ticks")
        p = sns.barplot(data=cog_descriptions, x="proportion", y="description", hue="cog")
        p.tick_params(labelsize=16)
        p.legend_.remove()
        p.set_ylabel("")
        p.set_xlabel("Proportion of sequences in COG", fontsize=16)
        # Convert figure to image
        img_buf = BytesIO()
        plt.savefig(img_buf, dpi = 200, bbox_inches='tight')
        pdf.image(img_buf, w = 180, x = None, y = None, type = '', link = '')
        plt.close()
        pdf.write(5,"\n\n")

        ###### SUB-SUBTITLE -- DISTRIBUTION OF UPSTREAM REGION LENGTH ######
        pdf.subtitle2(f"{i}.2 Distribution of upstream region lengths in the OrthoGroup ({cog})")
        #### Description ####
        pdf.normal_text("The upstream region lengths were extracted from the genome sequences, based on the gene locations of genes in the same OrthoGroup as your target.")
        #### Plot Upstream Region Lengths ####
        ups_100, ups_300, ups_500 = get_upstream_lengths(tf_data, cog)
        ## plot the three histograms in one figure
        fig = plot_upstream_lengths(ups_100, ups_300, ups_500)
        # Convert figure to image
        img_buf = BytesIO()
        plt.savefig(img_buf, dpi = 200, bbox_inches='tight')
        pdf.image(img_buf, w = 180, x = None, y = None, type = '', link = '')
        plt.close()
        pdf.add_page()

        ###### SUB-SUBTITLE -- MEME results and MOTIF LOCATION DISTRIBUTION ######
        pdf.subtitle2(f"{i}.3 MEME results and Motif Location Distribution for ({cog})")
        #### Description ####
        pdf.normal_text("The following plots show the MEME results and the distribution of motif locations for the OrthoGroup.")
        #### Plot MEME Results and Motif Location Distribution ####
        cog_meme = cog2png_df[cog2png_df["cog"] == cog]
        if cog_meme.shape[0] == 0:
            pdf.normal_text(f"No MEME results found for COG {cog}, skipping to next target.")
            continue
        for index, row in cog_meme.iterrows():
            png_path1,png_path2,png_path3,png_path4 = row["png_path_1"], row["png_path_2"], row["png_path_3"], row["png_path_4"]
            tsv_path1,tsv_path2,tsv_path3,tsv_path4 = row["tsv_path_1"], row["tsv_path_2"], row["tsv_path_3"], row["tsv_path_4"]
            matrix_path1,matrix_path2,matrix_path3,matrix_path4 = row["matrix_path_1"], row["matrix_path_2"], row["matrix_path_3"], row["matrix_path_4"]
            ups = row["ups"]
            mode = row["mode"]
            mode_long = "Any Number of Repetitions" if mode == "anr" else "Zero or One Occurrence Per Sequence"
            motiflength = row["motiflength"]
            #### Sub-sub-subtitle per MEME result ####
            pdf.meme_section(f"MOTIFS - {ups}bp upstream, {mode.upper()}, max length {motiflength}")
            #### Make the motif location distribution plot ####
            tsvs = [pd.read_csv(tsv_path, sep = "\t", header = 0) for tsv_path in [tsv_path1, tsv_path2, tsv_path3, tsv_path4] if tsv_path is not None and os.path.exists(tsv_path)]
            if not tsvs:
                pdf.normal_text(f"No valid TSV files found for COG {cog}, skipping to next target.")
                continue
            tsv1, tsv2, tsv3, tsv4 = tsvs[0], tsvs[1], tsvs[2], tsvs[3]
            fig = plot_meme_locations(tsv1, tsv2, tsv3, tsv4, png_path1, png_path2, png_path3, png_path4, ups)
            img_buf = BytesIO()
            plt.savefig(img_buf, dpi = 200, bbox_inches='tight')
            pdf.image(img_buf, w = 180, x = None, y = None, type = '', link = '')
            plt.close()
            pdf.write(5,"\n\n")
            pdf.add_page()
        pdf.write(5,"\n\n")
        i += 1

    pdf.output(f"{outdir}/target_protein_report.pdf")


def parsing_arguments():
    parser = argparse.ArgumentParser(description="Generate a report of target proteins from orthogroup fasta files.")
    parser.add_argument("-t", "--target_proteins", required=True, help="Comma-separated list of target protein IDs.")
    parser.add_argument("--tf_data", required=True, help="Path to the tf_data directory, created during OphiAURTHO pipeline.")
    parser.add_argument("--meme_info", required=True, help="Path to the meme_info directory, created during OphiAURTHO pipeline.")
    parser.add_argument("-o", "--outdir", required=True, help="Path to output directory for the reports.")

    args = vars(parser.parse_args())
    targets = args["target_proteins"].split(",")
    tf_data = os.path.abspath(args["tf_data"])
    meme_info = os.path.abspath(args["meme_info"])
    outdir = os.path.abspath(args["outdir"])

    return targets, tf_data, meme_info, outdir

def parse_location(location):
    """
    parse_location parses the output to the regex match on a cds.faa header, to extract the location information. 
    This can look like: 
    On positive strand: normal [11178..11507] || join(11854..12282,12340..12810)
    On negative strand: [complement(876..1205)] || complement(join(3015..3209,3255..3981,4032..5107))
    If the location contains "<" or ">", means it's a partial gene. Skipping those for now.
    """
    coords_str = re.findall(r'\d+', location) # returns string
    coords = list(map(int, coords_str))
    #print(coords) # [181763, 181846, 181897, 182256, 182311, 182691] -> group by 2 and then get range(2nd+1,3rd-1) etc
    if re.search(r"^\d+\.\.\d+$|^join\(\d+.*\)", location):
        strand = "+"
        start = min(coords)
        end = max(coords)
    if re.search(r"^complement\(\d+\.\.\d+\)|^complement\(join\(\d+.+\)", location):
        strand = "-"
        start = max(coords)
        end = min(coords)
    return start, end, strand

def get_cog_seqids(tf_data):
    """
    This function reads each Orthogroup.fasta file, and returns a list of seqids (gene_number) contained in each Orthogroup.
    """
    ## Find all COGs in the tf_data/cogs directory
    fastaList = glob.glob(f"{tf_data}/*/cogs/*.fasta")
    with open(f"{tf_data}/tf_cog2seqid_location.csv", mode = 'w') as out:
        out.write("cog;id;start;end;strand\n")
    for fasta in fastaList:
        cog = os.path.basename(fasta).replace(".fasta", "")
        with open(fasta, mode = 'rt') as f:
            fastdic = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        for key, value in fastdic.items():
            description = value.description
            location = re.search(r"location=(.*)\]{1}? ", description)[1]
            if re.search(r"<|>", location): continue
            start, end, strand = parse_location(location)
            with open(f"{tf_data}/tf_cog2seqid_location.csv", mode = 'a') as out:
                out.write(f"{cog};{key};{start};{end};{strand}\n")
    # read the file and return a pandas DataFrame
    cog_seqids = pd.read_csv(f"{tf_data}/tf_cog2seqid_location.csv", sep = ";", header = 0)
    cog_seqids['GeneNumber'] = cog_seqids["id"].apply(lambda x: x.split("|")[0])
    cog_seqids['ProtID'] = cog_seqids["id"].apply(lambda x: x.split("|")[1] if "|" in x else x)
    # Drop the start, end and strand columns, as they are not needed for the final output
    cog_seqids.drop(columns=["start", "end", "strand"], inplace=True)
    # write new file
    cog_seqids.columns = ["cog", "geneid", "GeneID", "ProtID"]
    cog_seqids.to_csv(f"{tf_data}/tf_cog2seqid.tsv", sep = "\t", index = False)
    return cog_seqids

def get_protein_descriptions(cog, tf_data):
    """
    This function reads a given Orthogroup.fasta file, and returns a dataframe with the protein descriptions and their counts.
    """
    tf_df = pd.DataFrame()
    tf_fam = cog.split("_")[0]
    fasta = f"{tf_data}/{tf_fam}/cogs/{cog}.fasta"
    # open file and keep only lines starting with >
    headers = os.popen(f"grep '>' {fasta}").read().split("\n")
    # extract "[protein=...]" from headers
    headers = [h.split("[protein=")[1].split("]")[0] for h in headers if "[protein=" in h]
    # turn to table, couting occurences of each header
    headers = pd.Series(headers).value_counts().reset_index()
    headers.columns = ["description", "count"]
    # add a cog column in first position
    headers.insert(0, "cog", cog)
    headers.insert(1, "tf_fam", tf_fam)
    tf_df = pd.concat([tf_df, headers], ignore_index=True)

    tf_df['total'] = tf_df.groupby('cog')['count'].transform('sum')
    tf_df['proportion'] = tf_df['count'] / tf_df['total']
    
    return tf_df

def get_upstream_lengths(tf_data, cog):
    # define dataframes for each upstream length
    ups_100 = pd.DataFrame()
    ups_300 = pd.DataFrame()
    ups_500 = pd.DataFrame()

    tf_fam = cog.split("_")[0]
    ups100 = f"{tf_data}/{tf_fam}/upstream/{cog}_upstream_100.fasta"
    ups300 = f"{tf_data}/{tf_fam}/upstream/{cog}_upstream_300.fasta"
    ups500 = f"{tf_data}/{tf_fam}/upstream/{cog}_upstream_500.fasta"

    # read the fasta with biopython
    ups100_records = [[rec.id, len(rec.seq)] for rec in SeqIO.parse(ups100, "fasta")]
    ups300_records = [[rec.id, len(rec.seq)] for rec in SeqIO.parse(ups300, "fasta")]
    ups500_records = [[rec.id, len(rec.seq)] for rec in SeqIO.parse(ups500, "fasta")]
    # make into df
    ups100_df = pd.DataFrame(ups100_records, columns = ["id", "length"]).assign(tf_fam = tf_fam, cog = cog, extract_length = "100bp")
    ups300_df = pd.DataFrame(ups300_records, columns = ["id", "length"]).assign(tf_fam = tf_fam, cog = cog, extract_length = "300bp")
    ups500_df = pd.DataFrame(ups500_records, columns = ["id", "length"]).assign(tf_fam = tf_fam, cog = cog, extract_length = "500bp")
    # concatenate to main df
    ups_100 = pd.concat([ups_100, ups100_df], ignore_index=True)
    ups_300 = pd.concat([ups_300, ups300_df], ignore_index=True)
    ups_500 = pd.concat([ups_500, ups500_df], ignore_index=True)
    # return the three dataframes
    return ups_100, ups_300, ups_500

def get_cog2png_df(tf_data, meme_info, cogs):
    # make a dictionary for each cog --> meme params and ups length --> png path
    cog2png_df = pd.DataFrame()
    for cog in cogs:
        tf_fam = cog.split("_")[0]
        for ups in [100, 300, 500]:
            for mode in ["anr", "zoops"]:
                for motiflength in [10, 20, 30]:
                    for nmotifs in [1,2,3,4]:
                        png_path = os.path.join(tf_data, f"{tf_fam}/meme/meme_{cog}_{ups}_{mode}_{motiflength}/logo_rc{nmotifs}.png")
                        tsv_path = os.path.join(meme_info, f"{cog}_{ups}_{mode}_{motiflength}_motif{nmotifs}.tsv")
                        matrix_path = os.path.join(meme_info, f"{cog}_{ups}_{mode}_{motiflength}_motif{nmotifs}.fasta")

                        if os.path.exists(png_path) and os.path.exists(tsv_path):
                            cog2png_df = pd.concat([cog2png_df, pd.DataFrame([[cog, ups, mode, motiflength, nmotifs, png_path, tsv_path, matrix_path]], columns = ["cog", "ups", "mode", "motiflength", "nmotifs", "png_path", "tsv_path", "matrix_path"])], ignore_index=True)
    # Pivot the df longer, to have one row per cog,ups,mode,motiflength, and 1 column per tsv_path_nmotif, png_path_nmotif and matrix_path_nmotif
    # final columns: cog, ups, mode, motiflength, tsv_path_1, png_path_1, matrix_path_1, tsv_path_2, png_path_2, matrix_path_2, tsv_path_3, png_path_3, matrix_path_3, tsv_path_4, png_path_4, matrix_path_4
    cog2png_df = cog2png_df.pivot_table(index=["cog", "ups", "mode", "motiflength"], columns="nmotifs", values=["tsv_path", "png_path", "matrix_path"], aggfunc="first")
    cog2png_df.columns = [f"{col[0]}_{col[1]}" for col in cog2png_df.columns]
    cog2png_df = cog2png_df.reset_index()
    
    return cog2png_df

def plot_upstream_lengths(ups_100, ups_300, ups_500):
    fig, axs = plt.subplots(1, 3, figsize=(10, 5))
    sns.set_style("white")
    sns.set_style("ticks")
    g_100 = sns.histplot(ups_100, x="length", multiple="dodge", binwidth=10, binrange = (1,111), kde=True, ax=axs[0])
    g_100.set_title("100bp upstream", fontsize=16)
    g_100.set_xlabel("Length (bp)", fontsize=14)
    g_100.set_ylabel("Count", fontsize=14)
    g_100.tick_params(labelsize=14)
    g_300 = sns.histplot(ups_300, x="length", multiple="dodge", binwidth=10, binrange = (1,311), kde=True, ax=axs[1])
    g_300.set_title("300bp upstream", fontsize=16)
    g_300.set_xlabel("Length (bp)", fontsize=14)
    g_300.set_ylabel("", fontsize=14)
    g_300.tick_params(labelsize=14)
    g_500 = sns.histplot(ups_500, x="length", multiple="dodge", binwidth=10, binrange = (1,511), kde=True, ax=axs[2])
    g_500.set_title("500bp upstream", fontsize=16)
    g_500.set_xlabel("Length (bp)", fontsize=14)
    g_500.set_ylabel("", fontsize=14)
    g_500.tick_params(labelsize=14)
    plt.tight_layout()
    return fig

def plot_meme_locations(tsv1, tsv2, tsv3, tsv4, png_path1, png_path2, png_path3, png_path4, ups):
    fig, axs = plt.subplots(nrows = 4, ncols = 2, figsize=(8, 8))
    # Reverse x axis for all plots (so that 1 is on the right and ups is on the left) and set it to negative values
    plot1, plot2, plot3, plot4 = sns.histplot(tsv1, x="Start", multiple="dodge", binwidth=10, binrange = (1,ups+1), kde=True, ax=axs[0,0]), sns.histplot(tsv2, x="Start", multiple="dodge", binwidth=10, binrange = (1,ups+1), kde=True, ax=axs[0,1]), sns.histplot(tsv3, x="Start", multiple="dodge", binwidth=10, binrange = (1,ups+1), kde=True, ax=axs[2,0]), sns.histplot(tsv4, x="Start", multiple="dodge", binwidth=10, binrange = (1,ups+1), kde=True, ax=axs[2,1])
    plot1.set_title(f"Motif-1", fontsize=16)
    plot1.set_xlabel("Position in upstream region (bp)", fontsize=14)
    plot1.set_ylabel("Count", fontsize=14)
    plot1.tick_params(labelsize=14)
    plot2.set_title(f"Motif-2", fontsize=16)
    plot2.set_xlabel("Position in upstream region (bp)", fontsize=14)
    plot2.set_ylabel("", fontsize=14)
    plot2.tick_params(labelsize=14)
    plot3.set_title(f"Motif-3", fontsize=16)
    plot3.set_xlabel("Position in upstream region (bp)", fontsize=14)
    plot3.set_ylabel("Count", fontsize=14)
    plot3.tick_params(labelsize=14)
    plot4.set_title(f"Motif-4", fontsize=16)
    plot4.set_xlabel("Position in upstream region (bp)", fontsize=14)
    plot4.set_ylabel("", fontsize=14)
    plot4.tick_params(labelsize=14)
    img1, img2, img3, img4 = mpimg.imread(png_path1), mpimg.imread(png_path2), mpimg.imread(png_path3), mpimg.imread(png_path4)
    axs[1,0].imshow(img1)
    axs[1,0].axis('off')  # turn off axis
    axs[1,1].imshow(img2)
    axs[1,1].axis('off')  # turn off axis
    axs[3,0].imshow(img3)
    axs[3,0].axis('off')  # turn off axis
    axs[3,1].imshow(img4)
    axs[3,1].axis('off')  # turn off axis
    plt.tight_layout()

    return fig

if __name__ == "__main__":
    main()