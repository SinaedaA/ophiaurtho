
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils
import gzip
import re
import argparse
import glob
import os

def main():
    fna_dir, gff_dir, create_db, tf_dir, length, stop_cds, logfile = parsing_arguments()
    query_dir = os.path.abspath(f"{tf_dir}/cogs")
    ups_dir = os.path.abspath(f"{tf_dir}/upstream")
    if not os.path.exists(ups_dir):
        print(f"Creating upstream directory {ups_dir}...")
        os.makedirs(ups_dir)
    if not os.path.exists(query_dir):
        print(f"Query directory {query_dir} does not exist, exiting...")
        exit(1)
    
    ## 1. Initiate logfile in 'w' mode to overwrite previous one
    with open(logfile, mode = 'w') as log:
        log.write("Program for extracting upstream regions of given COG files (fasta).")
        log.write(f"Options entered: \n--fna-dir {fna_dir} --gff-dir {gff_dir} --create-db {create_db} --query-dir {query_dir} --length {length} --stop-cds {stop_cds} --logfile {logfile}")
    
    ## 2. Get the fna dictionary (keys = contigs, values = seqrecords)
    fna_dict = fna_dir_to_dict(fna_dir)
    
    ## 3. Get reduced gene and cds dictionary (keys = locus_tag / accession, values = contig, start, end, strand)
    reduced_gene_dict, reduced_cds_dict = gff_dir_to_dict(gff_dir, create_db)
    #print(reduced_gene_dict) # keys for both are the gene_number duo
    
    ## 4. Glob query files
    queries = glob.glob(f"{query_dir}/*.fasta")
    ## 5. For each query file
    for query in queries:
        ## 5.1. Create output directory if it doesn't exist
        cog = query.split("/").pop().split(".")[0]
        out = f"{ups_dir}/{cog}_upstream"
        
        ## 5.2. Create oufile name
        outfile = f"{cog}_upstream_{length}.fasta"
        
        ## 5.3. Create query dictionary (keys = accession, values = dict(contig, locus_tag, start, end, strand))
        query_dict = make_query_dict(query)
        
        ## 5.4. Extract upstream sequences and write to logfile
        upstream_seqs = get_upstream_sequence(query_dict, fna_dict, reduced_gene_dict, stop_cds, length, logfile)
        
        ## 5.5. Write to outfile (fasta)
        with open(f"{ups_dir}/{outfile}", mode = 'w') as out:
            SeqIO.write(upstream_seqs, out, 'fasta')
        

def parsing_arguments():
    parser = argparse.ArgumentParser(
    prog='fetch_ups.py',
    description='Fetches upstream or regulatory regions from query CDSs.')
    parser.add_argument('--fna-dir',
                        help="Path to the directory Species.fna files (genomes).")
    parser.add_argument('--gff-dir', 
                        help = "Path to the directory with the Species.gff(.gz) files (genome annotations). Will create db in same directory")
    parser.add_argument('--create-db', default = True,
                        help = "Whether to create Species.gff.db files. If set to False, and no db files exist in gff-dir, error will be raised.")
    parser.add_argument('--tf-dir', default = "./",
                        help="Path to the TF directory. Should contain a 'cogs' directory, created during TF classficiation, qith query files (fasta files containing genes for which to extract regulatory regions).")
    parser.add_argument('--length', default = '1000', 
                        help = "Length of the upstream region to extract.")
    parser.add_argument('--stop-cds', default = True,
                        help = "Whether to stop the sequence extraction if another CDS is encountered before the end of the upstream length.")
    parser.add_argument('--logfile', default = "./logfile", 
                        help="Log file to show if index.cds.faa has been created")

    args = vars(parser.parse_args())
    fna_dir = os.path.abspath(args["fna_dir"])
    gff_dir = os.path.abspath(args["gff_dir"])
    create_db = args["create_db"]
    tf_dir = os.path.abspath(args["tf_dir"])
    length = int(args["length"])
    stop_cds = args["stop_cds"]
    logfile = os.path.abspath(args["logfile"])

    return fna_dir, gff_dir, create_db, tf_dir, length, stop_cds, logfile

def fna_dir_to_dict(fna_dir):
    """
    This function takes a directory containing fasta files, then reads the content of each and stores it in a dictionary, which it returns.
    It can handle gz zipped files, as well as unzipped fasta files, and doesn't care about the used alphabet.
    The keys are the contig names, and the values are the corresponding SeqRecords.
    """
    fnas = glob.glob(f"{fna_dir}/*fna*")
    fna_dict = {}
    for fna in fnas:
        if fna.endswith(".gz"):
            with gzip.open(fna, mode='rt') as zipf:
                content = SeqIO.to_dict(SeqIO.parse(zipf, "fasta"))
        else:
            with open(fna, mode='rt') as f:
                content = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        fna_dict.update(content)
    return fna_dict

def gff_dir_to_dict(gff_dir, create_db):
    """
    This function takes a directory containing the strain.gff and/or strain.gff.gz. 
    If create_db == True, will create strain.gff.db files in gff_dir, which can be loaded in further runs.
    1. Create gff_dict and gff_cds_dict, to further create reduced dictionaries. 
        - gff_dict: Keys = strains, values = db (features = gene) (gffutils)
        - gff_cds_dict: Keys = strains, values = db (features = CDS) (gffutils)
    2. Reduced dictionaries for genes and CDSs are created and returned
        - reduced_gene_dict : keys = locus_tag, values = dict(contig, start, end, strand))
        - reduced_cds_dict : keys = accession, values = dict(contif, start, end)
    """
    gff_dict = {}
    gff_cds_dict = {}
    if create_db == 'True':
        print(f"Creating GFF databases in {gff_dir}. If dbs already exist, will overwrite them. You have 30 seconds to cancel if you want.")
        #time.sleep(30)
        ## can't glob.glob(*.gff.*) because if db has already been created, then it will also match *.gff.db, and that will raise error
        gffs = glob.glob(f"{gff_dir}/*.gff.gz") + \
            glob.glob(f"{gff_dir}/*.gff")
        for gff in gffs:
            path = gff.split(".")[0]
            key = gff.split("/")[-1].split(".")[0] # strain number
            db = gffutils.create_db(data=gff,
                                    #dbfn=":memory:",
                                    dbfn=f"{path}.gff.db",
                                    merge_strategy='create_unique',
                                    disable_infer_transcripts=True,
                                    disable_infer_genes=True,
                                    force = True,
                                    keep_order = True)
            gff_dict[key] = db.all_features(featuretype="gene")
            gff_cds_dict[key] = db.all_features(featuretype="CDS")
    else:
        print(f"Looking for GFF databases in {gff_dir}")
        gff_dbs = glob.glob(f"{gff_dir}/*.gff.db")
        for gff_db in gff_dbs:
            key = gff_db.split("/")[-1].split(".")[0]
            gff = gffutils.FeatureDB(gff_db, keep_order = True)
            gff_dict[key] = gff.all_features(featuretype="gene")
            gff_cds_dict[key] = gff.all_features(featuretype="CDS")
            #print(gff.all_features(featuretype = "CDS")) ## generator
    ## create a reduced_gene_dict as well, with locus_tag -> contig, start, end, strand
    reduced_gene_dict = {}
    for generator in gff_dict.values():
        for feat in generator:
            # print(dir(feat)) # to get the attributes of feat
            contig = feat.seqid
            locus_tag = feat.id.split("-")[1] # remove the gene- part
            start = feat.start
            end = feat.end
            strand = feat.strand
            reduced_gene_dict[locus_tag] = {"contig": contig, "start": start, "end": end, "strand": strand}
    ## create a reduced_cds_dict as well, with accession -> contig, start, end
    reduced_cds_dict = {}
    for generator in gff_cds_dict.values():
        for feat in generator:
            contig = feat.seqid
            accession = feat.id.split("-")[1] # remove the cds- part
            start = feat.start
            end = feat.end
            reduced_cds_dict[accession] = {"contig": contig, "start": start, "end": end}
    return reduced_gene_dict, reduced_cds_dict
    
def make_query_dict(query):
    """
    This function takes a directory containing fasta files, reads the headers and extracts useful information.
    This can be COG files, or co-regulated gene files. Header must contain information of cds_from_genomic file type.
    It extracts: the seq_accession, the contig, the start and stop position, and the strand, returns a dictionary (keys = accession).
    """
    query_dict = {}
    for line in open(query).readlines():
        if re.match(">", line):
            accession = line.split(" ")[0].split(">")[1]
            contig = re.sub(r"_cds_.*\n", "",line.split("|")[2])
            location = re.search(r"location=(.*)\]{1}? ", line)[1]
            locus_tag = re.search(r"locus_tag=(.*?)\] \[{1}?", line)[1]
            if re.search(r"<|>", location): continue
            start, end, strand = parse_location(location)
            query_dict[accession] = {"contig": contig, "locus_tag": locus_tag, "start": start, "end": end, "strand": strand}
    return query_dict

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

def get_upstream_sequence(query_dict, fna_dict, reduced_gene_dict, stop_cds, length, logfile):
    '''
    Takes each gene in the query, and extracts the upstream sequence of given length. 
    The 'reduced_gene_dict' is used to check if upstream gene (or downstream depending on strand) overlap with upstream sequence.
    It will cut the sequence at the start/stop of upstream gene if stop_cds == True.
    Returns a list of SeqRecords for each gene in the query. 
    '''
    upstream_sequences = [] ## will be seqrecord list
    for accession, features in query_dict.items():
        ## I'm using locus_tag, which I use in "check_overlap", to query the reduced_gene_dict for it's keys kind of
        key_name = accession.split("|")[0] # only keep the gene_number part (8_4051|WP_189845047.1)
        contig, locus_tag, start, end, strand = features.values()
        # extract_start ALWAYS needs to be smaller than extract_end
        if strand == '+':
            extract_start = start - length + 1
            extract_end = start - 1
            ## create upstream_region_range (ups_range) then check if overlap with previous gene
            ups_range = range(extract_start, extract_end)
            overlap, next_range = check_overlap(dic = reduced_gene_dict, lt = key_name, ups_range = ups_range)
            ## change extract_start if there is overlap and stop_cds is set to True in parameters
            if (overlap != None) and (stop_cds == 'True'):
                print("Overlap detected between the desired upstream region range, and the next or previous gene... \nChanging the extraction coordinates because --stop-cds = True")
                extract_start = overlap.start + 1
        if strand == '-':
            extract_start = start + 1
            extract_end = start + length + 1
            ## create upstream_region_range (ups_range) then check if overlap with previous gene
            ups_range = range(extract_start, extract_end)
            overlap, next_range = check_overlap(dic = reduced_gene_dict, lt = key_name, ups_range = ups_range)
            ## change extract_start if there is overlap and stop_cds is set to True in parameters
            if (overlap != None) and (stop_cds == 'True'):
                print("Overlap detected between the desired upstream region range, and the next or previous gene... \nChanging the extraction coordinates because --stop-cds = True")
                extract_end = overlap.start - 1
        upstream_sequence = fna_dict[contig].seq[extract_start:extract_end]
        record = SeqRecord(upstream_sequence, id = locus_tag, description = f"[accession: {accession}] [contig: {contig}] [length: {len(upstream_sequence)}bp] [coordinates: {extract_start}:{extract_end}]")
        ## Only append to upstream_sequences list, if the sequence is not of length 0 (cause MEME doesn't work if there are some headers that don't have a sequence)
        if len(upstream_sequence) > 0:
            upstream_sequences.append(record)
        else:
            continue
        ## printing to logfile
        with open(logfile, mode = 'a') as log:
            log.write(f"Upstream region of {locus_tag} has length {len(upstream_sequence)}: \nUpstream range: {ups_range}\nUpstream gene: {next_range}\nOverlap: {overlap}\nSeq:{upstream_sequence}")
    return upstream_sequences
        
def check_overlap(dic, lt, ups_range):
    '''
    Called by 'get_upstream_sequence'.
    Creates ki and ik dictionary, based on reduced_gene_dict.
    ki: keys = keys of reduced_gene_dict, values = index
    ik: keys = index, values = keys of reduced_gene_dict.
    1. Get the index of the locus_tag of interest (query) and look up the index on the next one.
    2. Get the key name of the next one, based on the index found in step 1
    3. Find the 'start' and 'end' of the next locus_tag and create range()
    4. Check for overlap with the upstream_sequence range()
    5. Return overlap (which is a range of the intersection OR None) and next_range
    '''
    ## define key-index and index-key dictionary
    ki = dict()
    ik = dict()
    for i, k in enumerate(dic):
        ki[k] = i
        ik[i] = k
    ## get index of locus_tag of interest
    index_of_lt = ki[lt]
    if dic[lt]['strand'] == '+':
        next_lt = ik[index_of_lt - 1]
    elif dic[lt]['strand'] == '-':
        next_lt = ik[index_of_lt + 1]
    ## start is already always min, and end is already max coordinate
    next_range = range(dic[next_lt]["start"], dic[next_lt]["end"])
    overlap = range_intersect(next_range, ups_range)
    return overlap, next_range

def range_intersect(r1, r2):
    return range(max(r1.start,r2.start), min(r1.stop,r2.stop)) or None

if __name__ == '__main__':
    main()
