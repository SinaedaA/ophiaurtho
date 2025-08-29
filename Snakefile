import os
import glob
import platform

system = platform.system()  # Returns 'Linux', 'Darwin' (macOS), or 'Windows'
print(f"Detected system: {system}")
if system == "Windows":
    raise RuntimeError("This workflow is not supported on Windows. Please use Linux or macOS.")

# ---------------------- #
# -- Global variables -- #
# ---------------------- # 
assembly_info_path = "results/genomes/assembly.info"
project_name = config["proteinortho_project"]
tf_families_list = config["tf_families"].split(",") if "tf_families" in config else []
ups_lengths = config["upstream_region_length"] if "upstream_region_length" in config else 300
meme_modes = config["meme_modes"]
meme_lengths = config["meme_lengths"]

# ------------------------- #
# -- Checkpoint Handling -- #
# ------------------------- #
def get_strains(wildcards):
    checkpoint_output = checkpoints.find_strains.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        return [line.strip() for line in f]

def cds_faa_inputs(wildcards):
    return expand("results/cds_faa/{strain}.cds.faa", strain=get_strains(wildcards))
def cds_gff_inputs(wildcards):
    return expand("results/cds_gff/{strain}.cds.gff", strain=get_strains(wildcards))

def ipr_validate_inputs(wildcards):
    return expand("results/interproscan/{strain}.tsv", strain=get_strains(wildcards))

def get_dbcan_output(wildcards):
    """Function to get final outputs after checkpoint completes"""
    # This will be called after find_strains checkpoint
    checkpoint_output = checkpoints.find_strains.get(**wildcards).output[0]
    # Read the strain list
    with open(checkpoint_output) as f:
        strains = [line.strip() for line in f]
    # Return the expanded output file paths, which will be evaluated in the rule all.
    return expand("results/dbcan/{strain}/.dbcan_done", strain=strains)

# --------------- #
# -- Main Rule -- #
# --------------- #
rule all:
    input:
        "results/.cds_faa_created",
        #get_dbcan_output,
        "results/.all_dbcan_done",
        # "results/proteinortho/.snakemake_validate",
        # "results/proteinortho/.snakemake_transform_validate",
        # "results/proteinortho/.snakemake_grab_proteins_validate",
        "results/.all_cogs_renamed",
        "results/.all_iprscan_done",
        "results/.classify_tfs_done",
        "results/.extract_upstream_regions_done",
        "results/.all_meme_done",
        "results/.dbcan_tf_integrated"

# -------------------------------- #
# ------- Download genomes ------- #
# -------------------------------- #
rule download_genomes:
    output:
        assembly_info_path
    conda:
        "workflow/envs/aurtho.yml"
    params:
        outdir="results/genomes",
        genus=config["genus"],
        level=config["assemblylevel"],
        overwrite=config["overwrite"],
        db=config["download_genomes_db"],
        repo=config["repository"],
        max_genomes=config["max_genomes"]
    shell:
        "python workflow/scripts/1_download_genomes.py "
        "--outdir {params.outdir} "
        "--genus {params.genus} "
        "--level-of-assembly {params.level} "
        "--overwrite {params.overwrite} "
        "--database {params.db} "
        "--repository {params.repo} "
        "--max-genomes {params.max_genomes}"

## I'm actually not using this... And it's not updating the config file. 
checkpoint find_strains:
    input:
        assembly_info_path,
        "results/.symlinks_created"
    output:
        "results/strain_list.txt"
    run:
        strain_paths = glob.glob("results/genomes/GCF_*/[0-9]*.cds.faa")
        strain_ids = [os.path.basename(x).split(".")[0] for x in strain_paths]
        config["strains"] = strain_ids  # Update config
        with open("results/strain_list.txt", "w") as out:
            for sid in strain_ids:
                out.write(sid + "\n")

# -------------------------------------------- #
# ------- cds.faa and cds.gff creation ------- #
# -------------------------------------------- #
rule make_cds_faa_gff:
    input:
        assembly_info_path
    output:
        touch("results/.cds_faa_created")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        genome_dir="results/genomes"
    shell:
        "python workflow/scripts/2_make_cds_faa_gff.py "
        "--infofile {input} "
        "--genome_dir {params.genome_dir} "
        "--logfile {output} "

rule make_symlinks:
    input:
        assembly_info=assembly_info_path,
        cds_faa_created="results/.cds_faa_created"
    output:
        touch("results/.symlinks_created")
    params:
        genome_dir="results/genomes"
    shell:
        "mkdir -p results/cds_faa results/cds_gff results/fna && "
        "bash workflow/scripts/make_symlinks.sh {input.assembly_info} {params.genome_dir}"

rule make_gff_db:
    input:
        assembly_info_path,
        "results/.symlinks_created"
    output:
        touch("results/.gff_db_created")
    conda:
        "workflow/envs/aurtho.yml"
    shell:
        "python workflow/scripts/make_gff_db.py --gff-dir results/cds_gff/"

# ------------------ #
# -- Interproscan -- #
# ------------------ #
if system == "Darwin":
    rule interproscan:
        threads: 4
        input:
            "results/.symlinks_created",
            faa="results/cds_faa/{strain}.cds.faa"
        output:
            tsv="results/interproscan/{strain}.tsv"
        container:
            "docker://interpro/interproscan:5.75-106.0"
        params:
            #cpus=config["ipr_cpu"] or 4,
            applications=config["ipr_appl"]
        shell:
            "mkdir -p results/interproscan && "
            "export JAVA_OPTS='-Xmx8G' && "  # Limit to 8GB RAM
            "/opt/interproscan/interproscan.sh --input {input.faa} -appl {params.applications} -f tsv -o {output.tsv} --cpu {threads} --disable-precalc"

    rule validate_iprscan:
        input:
            ipr_validate_inputs
        output:
            "results/.all_iprscan_done"
        shell:
            "touch {output}"

# --------------------- #
# ------- dbCAN ------- #
# --------------------- #
rule run_dbcan:
    input:
        "results/.symlinks_created",
        faa="results/cds_faa/{strain}.cds.faa",
        gff="results/cds_gff/{strain}.cds.gff"
    output:
        touch("results/dbcan/{strain}/.dbcan_done")  # Validation file per strain
    conda:
        "workflow/envs/dbcan.yml"
    shell:
        "run_dbcan easy_substrate --mode protein "
        "--output_dir results/dbcan/{wildcards.strain} "
        "--input_raw_data {input.faa} "
        "--threads 8 "
        "--db_dir resources/dbcan/db/ "
        "--gff_type NCBI_prok "
        "--input_gff {input.gff} "

rule validate_dbcan:
    input:
        get_dbcan_output
    output:
        "results/.all_dbcan_done"
    shell:
        "touch {output}"

# ---------------------------- #
# ------- Proteinortho ------- #
# ---------------------------- #
rule run_proteinortho:
    input:
        "results/.symlinks_created",
        cds_faa=cds_faa_inputs  # Use function instead of glob
    output:
        main_tsv=f"results/proteinortho/{project_name}.proteinortho.tsv",
        info=f"results/proteinortho/{project_name}.info",
        blast_graph=f"results/proteinortho/{project_name}.blast-graph",
        graph=f"results/proteinortho/{project_name}.proteinortho-graph",
        validate="results/proteinortho/.snakemake_validate"
    conda:
        "workflow/envs/proteinortho.yml"
    params:
        cpus=config["proteinortho_cpu"]
    shell:
        "mkdir -p results/proteinortho && "
        "proteinortho6.pl -project={project_name} -singles -cpus={params.cpus} -selfblast -verbose {input.cds_faa} && touch {output.validate}; "
        "mv {project_name}.proteinortho.tsv {output.main_tsv}; mv {project_name}.proteinortho-graph {output.graph}; mv {project_name}.info {output.info}; mv {project_name}.blast-graph {output.blast_graph}"

rule transform_proteinortho:
    input:
        cds_faa=cds_faa_inputs,  # Use function instead of glob
        main_tsv=f"results/proteinortho/{project_name}.proteinortho.tsv",
        graph=f"results/proteinortho/{project_name}.proteinortho-graph"
    output:
        graph_summary=f"results/proteinortho/{project_name}.proteinortho-graph.summary",
        html=f"results/proteinortho/{project_name}.proteinortho.html",
        xml=f"results/proteinortho/{project_name}.proteinortho.xml",
        validate="results/proteinortho/.snakemake_transform_validate"
    conda:
        "workflow/envs/proteinortho.yml"
    shell:
        "proteinortho_summary.pl '{input.graph}' > '{output.graph_summary}' && "
        "proteinortho2html.pl {input.main_tsv} {input.cds_faa} > {output.html} && "
        "proteinortho2xml.pl {input.main_tsv} > {output.xml} && "
        "touch {output.validate}"

checkpoint grab_proteins:
    input:
        cds_faa=cds_faa_inputs,
        main_tsv=f"results/proteinortho/{project_name}.proteinortho.tsv"
    output:
        validate="results/proteinortho/.snakemake_grab_proteins_validate"
    params:
        cpus=config["proteinortho_cpu"],
        cog_dir=os.path.join("results/proteinortho", config["cog_directory"]),
        min_prot=config["minimum_cog_size"]
    conda:
        "workflow/envs/proteinortho.yml"
    shell:
        "mkdir -p {params.cog_dir} && "
        "proteinortho_grab_proteins.pl -tofiles={params.cog_dir} "
        "-cpus={params.cpus} -minprot={params.min_prot} "
        "{input.main_tsv} {input.cds_faa} && "
        "touch {output.validate}"

def get_orthogroups(wildcards):
    """Dynamically get list of groups after grab_proteins completes"""
    return glob_wildcards("results/proteinortho/cog/{project_name}.proteinortho.tsv.OrthoGroup{group}.fasta")[1]

rule rename_orthogroup_fastas:
    input:
        # define the input, by using the config value for the project name, as well as the group wildcard
        src=lambda wildcards: f"results/proteinortho/cog/{config['proteinortho_project']}.proteinortho.tsv.OrthoGroup{wildcards.group}.fasta"
    output:
        dest="results/proteinortho/cog/OrthoGroup{group}.fasta"
        #validate=temp("results/proteinortho/cog/.Orthogroup{group}.renamed.done")
    shell:
        "mv {input.src} {output.dest}"

rule validate_cog_renaming:
    input:
        "results/proteinortho/.snakemake_validate",
        "results/proteinortho/.snakemake_transform_validate",
        "results/proteinortho/.snakemake_grab_proteins_validate",
        expand("results/proteinortho/cog/OrthoGroup{group}.fasta", group=get_orthogroups)
    output:
        "results/.all_cogs_renamed"
    shell:
        "touch {output}"

# ---------------------------- #
# ------- Classify TFs ------- #
# ---------------------------- #
# rule classify_tfs:
#     input:
#         # should only be executed after dbcan AND grab_proteins are done
#         dbcan_done=get_dbcan_output,
#         cogs_renamed="results/proteinortho/cog/.all_renamed.done"
#     output:
#         touch("results/tf_data/.classify_tfs_done")
#     conda:
#         "workflow/envs/aurtho.yml"
#     params:
#         cog_dir=os.path.join("results/proteinortho", config["cog_directory"]),
#         perc_threshold=config["tf_perc_threshold"],
#         tf_families=config["tf_families"]
#     shell:
#         "python workflow/scripts/5_classify_tfs.py --dbcan_dir results/dbcan/ --cog_dir {params.cog_dir} --outdir results/tf_data/ --tf_families {params.tf_families} --perc_threshold {params.perc_threshold}"

rule classify_ipr_tfs:
    input:
        ipr_validate="results/.all_iprscan_done",
        cogs_renamed="results/.all_cogs_renamed"
    output:
        validate=touch("results/.classify_tfs_done")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        ipr_dir="results/interproscan/",
        cog_dir=os.path.join("results/proteinortho", config["cog_directory"]),
        perc_threshold=config["tf_perc_threshold"],
        tf_families=config["tf_families"]
    shell:
        "python workflow/scripts/5_classify_tfs_IPR.py --ipr_dir {params.ipr_dir} --tf_domains resources/tf_ipr_domains.tsv --cog_dir {params.cog_dir} --outdir results/tf_data/ --tf_families {params.tf_families} --perc_threshold {params.perc_threshold} && "
        "touch {output.validate}"

rule extract_upstream_regions:
    input:
        "results/.gff_db_created",
        tf_validate="results/.classify_tfs_done"
    output:
        validate=touch("results/tf_data/{tf}/.ups_done_{length}"), # when I use temp(), I need to ensure that the file is created in the shell command with touch
        logs=os.path.join("results/tf_data/{tf}/upstream_{length}.log")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        stop_cds=config["upstream_stop_cds"],
        length=lambda wildcards: wildcards.length,
        tf_dir=lambda wildcards: f"results/tf_data/{wildcards.tf}"
    shell:
        "python workflow/scripts/6_fetch_ups.py "
        "--fna-dir results/fna/ "
        "--gff-dir results/cds_gff/ "
        "--create-db False "
        "--tf-dir {params.tf_dir} "
        "--length {params.length} "
        "--stop-cds {params.stop_cds} "
        "--logfile {output.logs}" # && "
        #"touch {output.validate}"

rule validate_extract_ups:
    input:
        expand("results/tf_data/{tf}/.ups_done_{length}", tf=tf_families_list, length=ups_lengths)
    output:
        "results/.extract_upstream_regions_done"
    shell:
        "touch {output}"

# ------------------------ #
# ------- run MEME ------- #
# ------------------------ #
rule run_meme:
    input:
        ups_validate="results/.extract_upstream_regions_done"
        #tf_ups_files=lambda wildcards: f"results/tf_data/{wildcards.tf}/upstream/{wildcards.meme_target}_upstream_{wildcards.length}.fasta"
    output:
        validate=temp("results/tf_data/{tf}/meme/.meme_done_{motif_length}_{meme_mode}") # when I use temp(), I need to ensure that the file is created in the shell command with touch
    params:
        motif_length=lambda wildcards: wildcards.motif_length,
        meme_mode=lambda wildcards: wildcards.meme_mode,
        tf_dir=lambda wildcards: f"results/tf_data/{wildcards.tf}"
    shell:
        "bash workflow/scripts/meme_loop.sh {params.tf_dir}/upstream/ {params.meme_mode} {params.motif_length} {params.tf_dir} && "
        "touch {output.validate}"

rule validate_meme:
    input:
        expand("results/tf_data/{tf}/meme/.meme_done_{motif_length}_{meme_mode}",
               tf=tf_families_list,
               motif_length=meme_lengths,
               meme_mode=meme_modes)
    output:
        "results/.all_meme_done"
    shell:
        "touch {output}"

# ------------------------------------- #
# ------ Integrate dbCAN and TFs ------ #
# ------------------------------------- #
rule integrate_dbcan_tf:
    input:
        dbcan_done=get_dbcan_output,
        tf_validate="results/.classify_tfs_done"
    output:
        outfile=touch("results/.dbcan_tf_integrated")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        method=config["tf_dbcan_method"],
        extension=config["tf_dbcan_extension"],
        neighbours=config["tf_dbcan_neighbours"]
    shell:
        "python workflow/scripts/integrate_dbcan_tf.py " 
        "--dbcan_dir results/dbcan/ "
        "--gff_db_dir results/cds_gff/ "
        "--db_dir resources/dbcan/db/ "
        "--tf_classification results/tf_data/TF_classification.tsv "
        "--method {params.method} "
        "--extension {params.extension} "
        "--neighbours {params.neighbours} "
        "--outdir results/dbcan_tf/"

# ------------------------ #
# ------- Clean up ------- #
# ------------------------ #
# rule cleanup_validation_files:
#     input:
#         "results/tf_data/.extract_upstream_regions_done"
#     output:
#         "results/tf_data/.cleanup_done"
#     shell:
#         """
#         rm -f results/tf_data/*/.ups_done_*
#         touch {output}
#         """