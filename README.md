## Introduction

## 1. Dependencies
### 1.1. Mamba or similar python environment manager
I personally work with mamba, because environment resolution is faster, and errors are clearer when they appear. But this should work with micromamba, conda, miniconda, etc. 

### 1.2. Installing Snakemake
The only environment and package you need to install is snakemake itself. The other environment requirements for the different steps of the pipeline are taken care of by snakemake. 
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

### 1.3. Clone this repository locally
You can do this using `git` or `gh`. 

```bash
## move to the directory where you want to clone it
cd ~/Documents/snakemake
## using git
git clone https://github.com/SinaedaA/ophiaurtho.git
## using gh
gh repo clone SinaedaA/ophiaurtho
## go into that directory
cd ophiaurtho
```
Both of these commands should have the same effect: create a directory called ophiaurtho, containing everything that is necessary to run the pipeline. 

### 1.4. MacOS specific dependencies (brew)
For the Interproscan database download, we need md5sum to check the accuracy of the downloaded data, and this is not basally installed on MacOS. 
Additionally, during the pipeline I use the `rename` command, which has to be installed.

Both of these commands can be installed with Homebrew !

```bash
brew install coreutils # includes the md5sum command
brew install rename
```

### 1.5. Download databases
The pipeline depends on **dbCAN** and **InterProScan**, for which we need to download external databases. The first step is thus to create those and get them locally on your machine. There is a snakemake rule for this as well: 

```bash
snakemake -s workflow/rules/databases.smk --cores 8 --use-conda --rerun-incomplete
```

**Note**: if your machine runs on Linux, this rule will also install a local InterProScan version. 

### 1.6. InterProScan installation (MacOS)
InterProScan is only available for Linux machines (and is installed on Linux during in 1.4.). However, using a docker image, we can run it on Darwin as well. Therefore, for portability, InterProScan is always run through the docker image. Full preparation for the InterProScan docker install are found here: ./InterProScan_MacOS.md. 

## 2. Simple snakemake run
The simplest way to run ophiaurtho, is like this:
```bash
## Linux
snakemake --configfile config/config.yml --config max_genomes=10 --cores 8 --use-conda

## macOS (add apptainer arguments, to run InterProScan)
snakemake --configfile config/config.yml --config max_genomes=10 --cores 8 --use-conda --use-apptainer --apptainer-args "--bind $PWD/resources/interproscan/db/interproscan-5.75-106.0/data:/opt/interproscan/data"
```

Here's a breakdown of each argument.
- `--configfile` simply indicates to snakemake where to find the configuration file for the pipeline. This is where you can modify the parameters to the pipeline (see 4. Modifying the config file).
- `--config max_genomes=10` indicates to snakemake that it should disregard the max number of downloaded genomes specified in the config file, and take 10 instead. Note: using only 10 genomes of the same genus in this type of comparative genomics approach does limit its predictive power. The more genomes are used, the better the predictions, but also, the longer the pipeline will take to run. 
- `--cores 8` indicates that the pipeline can parallelize 8 jobs, if each job uses only 1 thread (more on this in 5. Optimizing run time).
- `--use-conda` tells snakemake to use the specified conda environments (specified in each rule, when applicable). 
- `--use-apptainer` : similar, but for the docker image of Interproscan (only on macOS). 
- `--apptainer-args` : specifies additional arguments for apptainer, in this case, how to name the interproscan database inside the image (only on macOS). 

You can try this run. Note that Interproscan and dbCAN (but especially the former) take quite a while to run. If you have a lot of computing power, you can optimize the parallelization and the number of CPUs given to Interproscan, to speed up the run. It is by default kept to a conservative amount, in order to be able to run it on a laptop. 

## 3. Output
The pipeline should create the following directory tree (tree results -L 1):
```
results
├── cds_faa
├── cds_gff
├── dbcan
├── dbcan_tf
├── fna
├── genomes
├── interproscan
├── proteinortho
├── resources
└── tf_data
```
- `results/genomes` : contains the downloaded genomes, as well as the info file (summary of downloaded genomes).
- `results/cds_faa` , `results/cds_gff` , `results/fna` : contain the cds.faa, cds.gff and fna links (to the actual files, located in results>genomes>sub-directories), for easy access. These are renamed based on a numeric identifier for each strain. 
results/proteinortho : contains the proteinortho run results, as well as a sub-directory called "cogs", which contain the OrthoGroup fasta files. 
results/dbcan , results/interproscan : contain the dbCAN and InterProScan run results, respectively.
- `results/tf_data` : contains the transcription factor(TF) classification and statistics files. TF classification is done based on the InterProScan prediction (if a certain DNA-binding domain is present, then the protein is classified as the associated TF family). Further detail regarding how TFs are classified can be found in 6. TF classification strategy. Additionally, it contains a TF directory for each family of interest. In the example config file, we set the "LacI" and "GntR" families. Each of these sub-directories contain 3 folders: cogs (with corresponding OrthoGroups), upstream (with upstream region fasta files) and meme  (with MEME results).
- `results/dbcan_tf` : contains the dbCAN data integrated with the TF data. For each "cluster of CAZymes" identified in the genomes, we now added the information of whether or not a TF is present in the cluster, based on our own TF classification.

For a more detailed view of the sub-directories, run this command tree results -L 2 .

## 4. Modifying the config file - Personalizing your run
The default `config.yml` file is located in the config directory. Here is a breakdown of the parameters for each rule. 
### 4.1. Downloading genomes
Genomes are downloaded from NCBI, using the ftp server. It takes in the following parameters:

```yaml
### Configuration for downloading genomes
genus: "Streptomyces"
assemblylevel: "complete" # Options: {complete,chromosome,contig,scaffold,all}
overwrite: "True" # Options: {True,False}
download_genomes_db: "refseq" # Options: {refseq,genbank}
repository: "bacteria" # Options: {bacteria,fungi,invertebrate,metagenomes,other,plant,protozoa,vertebrate_mammalian,vertebrate_other,viral}
max_genomes: "300"
```
_**Breakdown of the parameters:**_
- `genus` : the default genus is "Streptomyces", but you can change it to any genus. To filter for genomes of this genus, for now it is simply a text search of the "assembly_summary" file downloaded from NCBI. 
- `assemblylevel` : reflects the quality of the genome assembly. We suggest to leave this on "complete", except if you are working with a relatively unstudied genus, of which few Complete Genome Assemblies have been made. 
- `overwrite` : whether to overwrite previous assembly downloads. This applies to sub-directories inside the output dir (results/genomes). 
- `download_genomes_db` : which NCBI database to download genomes from. Default is "refseq", which we recommend (curated, non-redundant sequences, of high quality). 
- `repository` : which NCBI repository to download genomes from. Set to "bacteria" by default, as the pipeline was originally developed, and likely most appropriate for prokaryotic TFBS prediction. 

### 4.2. InterProScan
InterProScan is used to predict protein domains (through IPR identifiers). If you want other applications than Pfam and the default IPR domain prediction to be run, you can add them. However, this increases the runtime significantly, and for TF classification, we only use the IPR domains. 

```yaml
### Configuration for InterproScan
ipr_appl: "Pfam" #(in addition to the default InterproScan)
#ipr_cpu: 4 # Number of CPUs to use for InterproScan
```
- `ipr_appl` : which InterProScan applications to run. Default, only Pfam , in addition to the default InterProScan, which returns the IPR identifier for protein domains. 
- `ipr_cpu` : commented out. By default, this is 4 (internal parameter to interproscan). Can be changed for optimization, but then other parameters need to be changed ! See X. Balancing parallelization, cores, cpus, threads. 

### 4.3. Proteinortho
ProteinOrtho computes the clusters of orthologous proteins (OrthoGroups). Here are the modifiable options:

```yaml
### Configuration for proteinortho
proteinortho_project: "streptomyces"
proteinortho_cpu: 8
cog_directory: "cog"
minimum_cog_size: 4 # min number of proteins to be considered to create the COG fasta files
```

- `proteinortho_project` : changes the name of the prefix of proteinortho output files. 
- `proteinortho_cpu` : changes how many CPUs are given to proteinortho. Set to 8 by default. 
- `cog_directory` : name of the directory where OrthoGroup fasta files will be stored. This directory is created inside the results/proteinortho directory. 
- `minimum_cog_size` : min number of proteins present inside an OrthoGroup to write the corresponding OrthoGroup.fasta file. Set to 4. Setting it to <2 does not make sense, as with comparative genomics, we need at least 2 sequences for a comparison to be possible. In our proof-of-concept study, we found that 4 was the minimum number of upstream sequences to give to MEME, where we could find a reliable predicted TFBS.

### 4.4. Transcription factor family identification
Transcription factors are identified based on the InterProScan results. There is a hand-curated TF domain file, located in the resources directory (also found at the end of this README), which links each DNA-binding domain to a specific family of transcription factors. 

_**Note**_: Regardless of the below parameters, all TF families are identified and present in the "TF_classification.tsv" output file. 

These parameters only affect which families will be considered for further analysis (writing the OrthoGroup.fasta files, extracting the upstream regions, and performing MEME on these DNA sequences). 

```yaml
### TF families of interest
tf_families: "LacI,GntR"
tf_perc_threshold: 70 # Percentage of proteins in the COG, for the COG to be considered in following analyses
```
- `tf_families` : is a comma-separated list of TF family names. These are the available options: {AraC, ArgR, ArsR, AsnC, BirA, Crp, DeoR, Fis, Fur, GntR, IclR, LacI, LexA, LuxR, LysR, MarR, PadR, PhoB, SinR, TetR, TrmB}.
- `tf_perc_threshold` : is the minimum percentage of proteins in an OrthoGroup identified as this specific family, for the OrthoGroup to be considered in downstream analysis. Default is 70 percent. 

_Example_: OrthoGroup1234 contains 10 proteins. 50% were identified as LacI (5/10). This group would not be considered for further analysis. 

### 4.5. Upstream region extraction
The length of the upstream gene region will affect MEME results. In our proof-of-concept, we tested several UPS lengths, and the most prolific region in terms of TFBS identification was 300 nucleotides upstream of the start codon. However, this was specific to the LacI family, and different TFs might have different TFBS positionings.

```yaml
### Upstream region extraction
upstream_region_length: [100, 300, 500] # Length of the upstream region to extract
upstream_stop_cds: True # Options: {True,False} - If True, stop codons will be included in the upstream region
```
- `upstream_region_length` : how many nucleotides to extract, upstream of the TF gene coordinates. By default, we extract 3 different lengths, and predict motifs on all sets of sequences, for completeness. 
- `upstream_stop_cds` : whether or not to stop extracting nucleotides when encountering the beginning or the end of another gene. Default: True. If neighboring genes are very conserved in different strains, we might over-predict DNA motifs, which actually represent a conserved gene region. 

### 4.6. MEME parameters
MEME is used to identify over-represented DNA motifs in the upstream regions of genes belonging to the same OrthoGroup (in this case, transcription factor OrthoGroups). 
```yaml
### MEME parameters
meme_modes: ["anr", "zoops"]
meme_lengths: [30, 20, 10]
```

- `meme_modes` : by default the pipeline will run both "anr" and "zoops". This parameter defines how MEME identifies a motif. ANR stands for "Any Number of Repeats", meaning the motif can be found any number of times in the same upstream region. ZOOPS stands for "Zero Or One Per Sequence", meaning the motif is allowed to be found only once (or zero times) per upstream region. Sometimes, TFs use repeated motifs (where often one is less conserved/perfect than the other) in order to increase . 
- `meme_lengths` : length of the MEME motifs found. By default, we run MEME with 3 different lengths. If it is known for your TF family of interest the typical length of the motif, you can restrict this parameter. 


### 4.7. TF data and dbCAN integration
This script aims at obtaining additional information about the identified TFs, more specifically, if they are found inside gene clusters of CAZymes. To do this, we rely on the CGC prediction of dbCAN, and identify whether or not any of our TFs are located inside, or close to the cluster. There are 2 ways to add a TF to a CGC; the neighbouring method, and the overlap method (or both ). See below for details. 
Additionally, it incorporates substrate prediction data to a final CGC summary table, at the gene level, and the cluster level. Indeed, dbCAN predicts the substrate of each CAZyme that it identifies, and for each CGC, it tries to predict a cluster-level substrate, when possible. These data are located in different tables, which are assembled into one single table through this script.
```yaml
### Integrate dbCAN and TF data
tf_dbcan_method: "neighbouring" # Options: {neighbouring,overlap,both} - Method to integrate dbCAN and TF data
tf_dbcan_extension: "0" # Extension length for the CGC location in base pairs, to get neighbouring TFs
tf_dbcan_neighbours: "2" # Number of CGC neighbouring genes, to get neighbouring TFs
```

- `tf_dbcan_method`: which method to use for TF-CGC integration. 
- `overlap`: add a TF to the CGC if the gene coordinates of the TF overlap with the cluster coordinates of the CGC.
- `neighbouring`: add a TF to the CGC is the TF gene is a (direct or indirect) neighbour of this cluster.
- `both` : use both methods.
- `tf_dbcan_extension` : number of nucleotides to extend the CGC cluster limits by, on each side. Default is 0 (no extension). When the method is "neighbouring", this parameter is internally set to 0 again. 
- `tf_dbcan_neighbours` : number of neighbouring genes to consider in checking if they are previously identified TFs. Default is 2, meaning it will extract the gene_ids of the 2 genes directly preceding the start coordinate of the CGC, as well as those of the 2 gene following the end coordinate of the CGC, and check if any of these 4 total genes are present in our identified TF table. If so, these TFs will be added to the CGC, and the limits of the latter will be extended.

_**Note**_: when TFs are added to a CGC, the start and end coordinates of the latter are updated in the final output table.

This script returns a summary CGC table, containing for each CGC, the CAZymes and other genes, the newly added TFs, and the predicted substrate (for the CGC, and for each individual CAZyme). 

## 5. Optimizing run time

As mentioned, InterProScan is quite resource intensive, even in its most conservative tool usage. We can give it more CPUs, so it can analyse more sequences simultaneously, however this means we have to change some inputs we give to Snakemake, and some Lima configuration parameters. Here's a breakdown of why: 

### Snakemake cores and threads
Snakemake manages core usage globally, and you can specify how many cores you allow it to use in the global snakemake call. Example: 
```bash
snakemake --configfile config/config.yml --config max_genomes=10 --cores 8 --use-conda --use-apptainer --apptainer-args "--bind $PWD/resources/interproscan-5.75-106.0/data:/opt/interproscan/data"
```
With the `--cores` flag, we tell it it can use 8 cores. 

Inside the Snakefile, for each rule of the pipeline, we can also specify how many threads the rule uses. Let's assume each core has 1 thread for computation. If a rule uses only 1 thread (default), that means that snakemake can run 8 instances of the same rule at the same time (as it has 8 cores). This is the case for non-resource intensive rules. For example, the rule that creates the "cds.faa" files uses a python script that is very quick and doesn't need a lot of computation power. Therefore, it will run this rule on 8 genomes at the same time, 1 thread each, 1 core each. 

_**Note**_: if the `threads` setting of a rule is higher than the `--cores` flag, then it defaults to the latter. 

### Lima VM CPU and memory
Remember, during the pipeline set-up, we added this to the lima config file:
```yaml
### For CPU usage
cpus: 8
memory: "16GiB"
disk: "100GiB"
```
The `memory` line says that the VM will have access to 16GiB of RAM (active memory), to store currently used data.

The `disk` line gives the VM 100GiB of storage capacity. If you generate more than 100GiB of data, an error will occur, but this shouldn't be the case for Interproscan.

Finally the `cpus` parameter is set to 8. That means that, out of the CPU cores available on your machine, only 8 are made accessible to the VM. 

### InterProScan - CPUs parameter
InterProScan itself is run using Lima (`apptainer`), the virtual machine. The software has its own `--cpu` flag, where you can define how many cores it can use. In this case, it determines how many processes it can run simultaneously, or, how many sequences it can analyze in parallel. 

Imagine you tell InterProScan it can use 32 CPU cores, but the VM only has 4 CPU cores available to run, the software will run into an error (probably an OOM error). 

### InterProScan - Java memory
InterProScan is written in Java, and when running it, you can configure how much RAM it is allowed to use ! It is important that the amount of active memory used by InterProScan is within the limits of what is available to the Lima VM (set to 16GiB above). This is done through the `-Xmx` flag. 

### Summary
- **Snakemake** manages CPU cores globally; each rule requests threads.
- **InterProScan**’s `--cpu` controls internal parallelism and memory usage.
- **Java** memory (`-Xmx`) must be set below Lima VM’s available RAM.
- **Lima** VM must be configured with enough CPUs and memory.
- Balance `--cores`, `threads`, and InterProScan `--cpu` to optimize speed and avoid OOM kills.

### Current settings
- Snakemake runs with 8 cores (`--cores 8` in the global snakemake call)
- Interproscan rule runs with 4 threads (`threads: 4`). Only 2 Interproscan analyses can be done simultaneously (2*4 = 8).
- Lima VM has 8 cores available (`cpus: 8` in the "~/.lima/apptainer/lima.yaml", lima config file), so we can run 2x4 core jobs.
- Lima VM has 16GiB RAM available (`memory: 16GiB` , also in the lima config file)
- Java memory is set to 8GiB (`-Xmx8G`) inside the rule, meaning that 2 Interproscans can run on the 16GiB VM. 
- InterProScan runs with 4 cpus (`--cpus {threads}`), set inside the rule as "cpus = threads for that rule". 

As you can see, my Lima VM parameters are double what I give to InterProScan (through java and internal parameters), so I can run 2 jobs in parallel.

### Increasing the parallelization
If you have more computation power available, you can of course increase parallelization, which will decrease run-time. 

To check how much cores you have available on your computer you can do it like this:
```bash
## MacOS
sysctl hw.physicalcpu # actual number of cores available
sysctl hw.logicalcpu # logical number of cores available (threads)

## Linux
nproc # number of physical cores
lscpu # look for "Thread(s) per core"
```

_**Note**_: imagine you have 16 cores, it is not really a good idea to give the entirety of those 16 cores to the VM, as your computer needs some power to work as well. I once crashed an old computer doing this. Maybe now, computers have safety measures for this, but I haven't tested it. 

The easiest (and safest) way to increase the parallelization and speed up the pipeline run, is to keep approximately the same proportion I set. 

#### Increasing internal InterProScan parallelization
First, we could double each of the parameters ! That way, we increase the Interproscan internal parallelization, and memory usage (more protein sequences run simultaneously).

|    Level     | Parameter | Value |    Where    |
| :----------: | :-------: | :---: | :---------: |
|  snakemake   |  --cores  |  16   | global call |
|     rule     |  threads  |   8   |  Snakefile  |
|     lima     |   cpus    |  16   |  lima.yaml  |
|     lima     |  memory   |  32   |  lima.yaml  |
|     java     |   -Xmx    |  16   |  Snakefile  |
| interproscan |   --cpu   |   8   |  Snakefile  |

#### Increasing the number of genomes analysed simultaneously by individual InterProScan processes
We could also double all of the parameters, **except** for the java and interproscan levels. This increases how many genomes (proteomes) snakemake will start analysing in individual InterProScan runs (4 instead of 2). Each InterProScan run however still uses the same number of threads and active memory. 

|    Level     | Parameter | Value |    Where    |
| :----------: | :-------: | :---: | :---------: |
|  snakemake   |  --cores  |  16   | global call |
|     rule     |  threads  |   8   |  Snakefile   |
|     lima     |   cpus    |  16   |  lima.yaml  |
|     lima     |  memory   |  32   |  lima.yaml  |
|     java     |   -Xmx    |   8   |  Snakefile   |
| interproscan |   --cpu   |   4   |  Snakefile   |

### Benchmarking increased parallelization
Testing was done on Linux, so I didn't use the lima parameters. I tested it on 100 genomes each time.

#### Doubling everything
|    Level     | Parameter | Value |    Where    |
| :----------: | :-------: | :---: | :---------: |
|  snakemake   |  --cores  |  16   | global call |
|     rule     |  threads  |   8   |  Snakefile   |
|     java     |   -Xmx    |  16   |  Snakefile   |
| interproscan |   --cpu   |   8   |  Snakefile   |

Time for full run: 

#### Quadrupling everything
|    Level     | Parameter | Value |    Where    |
| :----------: | :-------: | :---: | :---------: |
|  snakemake   |  --cores  |  32   | global call |
|     rule     |  threads  |  16   |  Snakefile   |
|     java     |   -Xmx    |  32   |  Snakefile   |
| interproscan |   --cpu   |  16   |  Snakefile   |

Time for full run: 

## 6. TF identification strategy
Based on the IPR domains contained in the Interproscan result tables, we classify proteins into different TF families. 

The program loads the IPR predictions and the TF family definition table (see below), and filters the former, like this: 
- Remove hits that don't contain an IPR domain
- Keep hits with `Status == T`. This means that the IPR hit is quite reliable and the domain has been curated. I refrained from using the `e-value` as a hard threshold, as the correct threshold for reliability depends on the database size, which makes it difficult to use a 'one threshold fits for all' approach.
- For each protein, keep the 3 best hits (based on e-value). 
- Keep only IPR ids that are present in the TF family definition table.

That gives us a dataframe of proteins that are potentially transcription factors. 

_Next, we want to know if the TFs belong to OrthoGroups, where most other proteins are predicted as the same TF_?
So, the TFs are matched to their respective OrthoGroups, and the content of each OrthoGroup is evaluated. Statistics per potential TF family COG are also outputted, detailing the percentage of proteins in the OrthoGroup, identified as belonging to the same TF family. The user can input a minimum percentage, as a threshold. OrthoGroups that don't meet this threshold are not considered for further analysis.
 
Here are the available TF families (use the Short_name in the config.yaml file).

| Family_name | Short_name |                   Family_description                   | InterPro  |
| :---------: | :--------: | :----------------------------------------------------: | :-------: |
|    ArgR     |    ArgR    |          arginine repressor C-terminal domain          | IPR020899 |
|     Fur     |    Fur     |                 Fur C-terminal domain                  | IPR043135 |
|  TetR/AcrR  |    TetR    | TetR/AcrR transcriptional repressor, C-terminal domain | IPR049444 |
|  TetR/AcrR  |    TetR    |                 DNA-binding HTH domain                 | IPR001647 |
|    LysR     |    LysR    |                         HTH_1                          | IPR000847 |
|    SinR     |    SinR    |           SinR repressor/SinI anti-repressor           | IPR001387 |
|  AraC/XylS  |    AraC    |      AraC/XylS family transcriptional regulators       | IPR050204 |
|    GntR     |    GntR    |       Bacterial regulatory proteins, gntR family       | IPR000524 |
|  GerE/LuxR  |    LuxR    |        Transcription regulator LuxR, C-terminal        | IPR000792 |
|    PhoB     |    PhoB    |           OmpR/PhoB-type DNA-binding domain            | IPR001867 |
|    MarR     |    MarR    |           MarR family, MarR-type HTH domain            | IPR000835 |
|  AsnC/Lrp   |    AsnC    |                  AsnC-type HTH domain                  | IPR000485 |
|     Fis     |    Fis     |            DNA binding HTH domain, Fis-type            | IPR002197 |
|  LacI/GalR  |    LacI    |                  LacI-type HTH domain                  | IPR000843 |
|    ArsR     |    ArsR    |            HTH ArsR-type DNA-binding domain            | IPR001845 |
|     Crp     |    Crp     |                  Crp-type HTH domain                   | IPR012318 |
|    IclR     |    IclR    |        Transcription regulator IclR, N-terminal        | IPR005471 |
|    BirA     |    BirA    |    Biotin operon repressor, helix-turn-helix domain    | IPR004409 |
|    LexA     |    LexA    |           LexA repressor, DNA-binding domain           | IPR006199 |
|    TrmB     |    TrmB    |        Transcription regulator TrmB, N-terminal        | IPR002831 |
|    PadR     |    PadR    |        Transcription regulator PadR, N-terminal        | IPR005149 |
|    DeoR     |    DeoR    |                  DeoR-type HTH domain                  | IPR001034 |

If any users want to add TFs, they can always contact me and contribute to this very simplified database. 

## 7. Writing ALL Orthogroups to fasta files
By default, in order to save space, the pipeline only writes the Orthogroups identified as TFs in the dataset. 

```bash
snakemake -s workflow/rules/get_all_cogs.smk --configfile config/config.yml --config timestamp=$TIMESTAMP --cores 8 --use-conda --rerun-incomplete
```

```bash
snakemake results/.all_iprscan_done --configfile config/config.yml --config max_genomes=10 --cores 8 --use-conda --use-apptainer --apptainer-args "--bind $PWD/resources/interproscan-5.75-106.0/data:/opt/interproscan/data" --benchmark-extended
```