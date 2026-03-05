import argparse
import os

def parsing_arguments():
    parser = argparse.ArgumentParser(
        prog='run_proteinortho',
        description='Runs proteinortho6 software. Assumes it is already installed and in PATH. If not, you can install it with conda.')

    parser.add_argument('--project_name', required=True,
                        help="Name of the project, which will be used as a suffix for proteinortho output files.")
    parser.add_argument('--cpus', default = "4",
                        help="How many CPUs to use for proteinortho. Default is 4, which will be quite slow.")
    parser.add_argument('--cdsfaa_files', required=True, 
                        help="Location of the cds.faa files. Example: /path/to/cds_faa/*.cds.faa")
    parser.add_argument('--min_prot', default=2,
                        help="What is the minimum number of proteins in an OrthoGroup to be written to a FASTA file. [Default=2].")
    parser.add_argument('--cog_dir', default = "COG",
                        help="Location of the COG directory. [Default=COG].")
    parser.add_argument('--logfile', default = "proteinortho",
                        help="Name of the log file (can contain path to as well). Default is proteinortho_{suffix}.log.")

    args = vars(parser.parse_args())
    project = args["project_name"]
    cpus = args["cpus"]
    faa = args["cdsfaa_files"]
    min_prot = args["min_prot"]
    cog_dir = args["cog_dir"]
    logfile = args["logfile"]

    return project, cpus, faa, min_prot, cog_dir, logfile

def main():
    project, cpus, faa, min_prot, cog_dir, logfile = parsing_arguments()
    # check if cog_dir exists, if not create it
    if not os.path.exists(cog_dir):
        os.makedirs(cog_dir)
    # Run proteinortho in the command line
    os.system(f"proteinortho6.pl -project={project} -singles -cpus={cpus} -selfblast -verbose {faa} 2&>1 | tee {logfile}.log")
    # Make a summary file for the proteinortho run
    os.system(f"proteinortho_summary.pl '{project}.proteinortho-graph' > '{project}.proteinortho-graph.summary' 2&>1 | tee {logfile}_summary.log")
    # Make the html and xml output files
    os.system(f"proteinortho2html.pl {project}.proteinortho.tsv {faa} > {project}.proteinortho.html 2&>1 | tee {logfile}_2html.log")
    os.system(f"proteinortho2xml.pl {project}.proteinortho.tsv > {project}.proteinortho.tsv.xml 2&>1 | tee {logfile}_2xml.log")
    # Grab the proteins and write them to FASTA files
    os.system(f"proteinortho_grab_proteins.pl -tofiles={cog_dir} -cpus={cpus} -minprot={min_prot} {project}.proteinortho.tsv {faa} 2&>1 | tee {logfile}_grab_proteins.log")

if __name__ == '__main__':
    main()