import os
import glob
import platform
from datetime import datetime

## System check
system = platform.system()  # Returns 'Linux', 'Darwin' (macOS), or 'Windows'
print(f"Detected system: {system}")
if system == "Windows":
    raise RuntimeError("This workflow is not supported on Windows. Please use Linux or macOS.")

## Determine timestamp and output directory
timestamp = config.get("timestamp")
outdir = f"results/{timestamp}"
print(f"Output directory: {outdir}")

# ---------------------- #
# -- Global variables -- #
# ---------------------- # 
assembly_info_path = f"{outdir}/genomes/assembly.info"
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
    return expand(f"{outdir}/cds_faa/{{strain}}.cds.faa", strain=get_strains(wildcards))
def cds_gff_inputs(wildcards):
    return expand(f"{outdir}/cds_gff/{{strain}}.cds.gff", strain=get_strains(wildcards))

def ipr_validate_inputs(wildcards):
    return expand(f"{outdir}/interproscan/{{strain}}.tsv", strain=get_strains(wildcards))

def get_dbcan_output(wildcards):
    """Function to get final outputs after checkpoint completes"""
    # This will be called after find_strains checkpoint
    checkpoint_output = checkpoints.find_strains.get(**wildcards).output[0]
    # Read the strain list
    with open(checkpoint_output) as f:
        strains = [line.strip() for line in f]
    # Return the expanded output file paths, which will be evaluated in the rule all.
    return expand(f"{outdir}/dbcan/{{strain}}/.dbcan_done", strain=strains)

# --------------- #
# -- Main Rule -- #
# --------------- #
rule all:
    input:
        f"{outdir}/.cds_faa_created",
        #get_dbcan_output,
        f"{outdir}/.all_dbcan_done",
        # f"{outdir}/proteinortho/.snakemake_validate",
        # f"{outdir}/proteinortho/.snakemake_transform_validate",
        # f"{outdir}/proteinortho/.snakemake_grab_proteins_validate",
        f"{outdir}/.all_cogs_renamed",
        f"{outdir}/.all_iprscan_done",
        f"{outdir}/.classify_tfs_done",
        f"{outdir}/.extract_upstream_regions_done",
        f"{outdir}/.all_meme_done",
        f"{outdir}/.dbcan_tf_integrated",
        f"{outdir}/.meme_info_extracted",
        f"{outdir}/benchmarks/combined_benchmarks.tsv"

# -------------------------------- #
# ------- Download genomes ------- #
# -------------------------------- #
rule download_genomes:
    output:
        assembly_info_path,
        f"{outdir}/benchmarks/download_genomes/benchmark.tsv"
    conda:
        "workflow/envs/aurtho.yml"
    params:
        outdir=f"{outdir}/genomes",
        genus=config["genus"],
        level=config["assemblylevel"],
        overwrite=config["overwrite"],
        db=config["download_genomes_db"],
        repo=config["repository"],
        max_genomes=config["max_genomes"],
        required=config["required_genomes"] if "required_genomes" in config else []
    benchmark:
        f"{outdir}/benchmarks/download_genomes/benchmark.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        "python workflow/scripts/1_download_genomes.py "
        "--outdir {params.outdir} "
        "--genus {params.genus} "
        "--level-of-assembly {params.level} "
        "--overwrite {params.overwrite} "
        "--database {params.db} "
        "--repository {params.repo} "
        "--max-genomes {params.max_genomes} "
        "--required-genomes {params.required}"

## I'm actually not using this... And it's not updating the config file. 
checkpoint find_strains:
    input:
        assembly_info_path,
        f"{outdir}/.symlinks_created"
    output:
        f"{outdir}/strain_list.txt"
    run:
        strain_paths = glob.glob(f"{outdir}/genomes/GCF_*/[0-9]*.cds.faa")
        #print(f"Strain paths: {strain_paths}")
        strain_ids = [os.path.basename(x).split(".")[0] for x in strain_paths]
        config["strains"] = strain_ids  # Update config
        with open(f"{outdir}/strain_list.txt", "w") as out:
            for sid in strain_ids:
                out.write(sid + "\n")

# -------------------------------------------- #
# ------- cds.faa and cds.gff creation ------- #
# -------------------------------------------- #
rule make_cds_faa_gff:
    input:
        assembly_info_path
    output:
        touch(f"{outdir}/.cds_faa_created")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        genome_dir=f"{outdir}/genomes"
    benchmark:
        f"{outdir}/benchmarks/make_cds_faa_gff/benchmark.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        "python workflow/scripts/2_make_cds_faa_gff.py "
        "--infofile {input} "
        "--genome_dir {params.genome_dir} "
        "--logfile {output} "

rule make_symlinks:
    input:
        assembly_info=assembly_info_path,
        cds_faa_created=f"{outdir}/.cds_faa_created"
    output:
        touch(f"{outdir}/.symlinks_created")
    params:
        genome_dir=f"{outdir}/genomes"
    benchmark:
        f"{outdir}/benchmarks/make_symlinks/benchmark.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        f"mkdir -p {outdir}/cds_faa {outdir}/cds_gff {outdir}/fna && "
        "bash workflow/scripts/make_symlinks.sh {input.assembly_info} {params.genome_dir}"

rule make_gff_db:
    input:
        assembly_info_path,
        f"{outdir}/.symlinks_created"
    output:
        touch(f"{outdir}/.gff_db_created")
    conda:
        "workflow/envs/aurtho.yml"
    benchmark:
        f"{outdir}/benchmarks/make_gff_db/benchmark.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        f"python workflow/scripts/make_gff_db.py --gff-dir {outdir}/cds_gff/"

# ------------------ #
# -- Interproscan -- #
# ------------------ #
if system == "Darwin":
    rule interproscan:
        threads: 4
        input:
            f"{outdir}/.symlinks_created",
            faa=f"{outdir}/cds_faa/{{strain}}.cds.faa"
        output:
            tsv=f"{outdir}/interproscan/{{strain}}.tsv"
        container:
            "docker://interpro/interproscan:5.75-106.0"
        params:
            #cpus=config["ipr_cpu"] or 4,
            applications=config["ipr_appl"]
        benchmark:
            f"{outdir}/benchmarks/interproscan/{{strain}}.tsv"
            #repeat("benchmarks/somecommand/{sample}.tsv", 3)
        shell:
            f"mkdir -p {outdir}/interproscan && "
            f"export JAVA_OPTS='-Xmx8G' && "  # Limit to 8GB RAM
            f"/opt/interproscan/interproscan.sh --input {input.faa} -appl {params.applications} -f tsv -o {output.tsv} --cpu {threads} --disable-precalc"

    rule validate_iprscan:
        input:
            ipr_validate_inputs
        output:
            validate=f"{outdir}/.all_iprscan_done"
        shell:
            "touch {output.validate}"

if system == "Linux":
    rule interproscan:
        threads: 4
        input:
            f"{outdir}/.symlinks_created",
            faa=f"{outdir}/cds_faa/{{strain}}.cds.faa"
        output:
            tsv=f"{outdir}/interproscan/{{strain}}.tsv"
        params:
            #cpus=config["ipr_cpu"] or 4,
            applications=config["ipr_appl"]
        benchmark:
            f"{outdir}/benchmarks/interproscan/{{strain}}.tsv"
            #repeat("benchmarks/somecommand/{sample}.tsv", 3)
        shell:
            f"mkdir -p {outdir}/interproscan && "
            f"export JAVA_OPTS='-Xmx8G' && "  # Limit to 8GB RAM
            f"./resources/interproscan/interproscan.sh --input {input.faa} -appl {params.applications} -f tsv -o {output.tsv} --cpu {threads} --disable-precalc"

    rule validate_iprscan:
        input:
            ipr_validate_inputs
        output:
            f"{outdir}/.all_iprscan_done"
        shell:
            "touch {output}"

# --------------------- #
# ------- dbCAN ------- #
# --------------------- #
rule run_dbcan:
    input:
        f"{outdir}/.symlinks_created",
        faa=f"{outdir}/cds_faa/{{strain}}.cds.faa",
        gff=f"{outdir}/cds_gff/{{strain}}.cds.gff"
    output:
        touch(f"{outdir}/dbcan/{{strain}}/.dbcan_done")  # Validation file per strain
    conda:
        "workflow/envs/dbcan.yml"
    benchmark:
        f"{outdir}/benchmarks/dbcan/{{strain}}.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        f"run_dbcan easy_substrate --mode protein "
        f"--output_dir {outdir}/dbcan/{{wildcards.strain}} "
        f"--input_raw_data {input.faa} "
        f"--threads 8 "
        f"--db_dir resources/dbcan/db/ "
        f"--gff_type NCBI_prok "
        f"--input_gff {input.gff} "

rule validate_dbcan:
    input:
        get_dbcan_output
    output:
        f"{outdir}/.all_dbcan_done"
        # if the interproscan concatenation of benchmark works, do it here as well
    shell:
        "touch {output}"

# ---------------------------- #
# ------- Proteinortho ------- #
# ---------------------------- #
rule run_proteinortho:
    input:
        f"{outdir}/.symlinks_created",
        cds_faa=cds_faa_inputs  # Use function instead of glob
    output:
        main_tsv=f"{outdir}/proteinortho/{project_name}.proteinortho.tsv",
        info=f"{outdir}/proteinortho/{project_name}.info",
        blast_graph=f"{outdir}/proteinortho/{project_name}.blast-graph",
        graph=f"{outdir}/proteinortho/{project_name}.proteinortho-graph",
        html=f"{outdir}/proteinortho/{project_name}.proteinortho.html",
        graph_summary=f"{outdir}/proteinortho/{project_name}.proteinortho-graph.summary",
        validate=f"{outdir}/proteinortho/.snakemake_validate"
    conda:
        "workflow/envs/proteinortho.yml"
    params:
        cpus=config["proteinortho_cpu"]
    benchmark:
        f"{outdir}/benchmarks/proteinortho/run_benchmark.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        f"mkdir -p {outdir}/proteinortho && "
        f"proteinortho6.pl -project={project_name} -singles -cpus={params.cpus} -selfblast -verbose {input.cds_faa} && touch {output.validate}; "
        f"mv {project_name}.proteinortho.tsv {output.main_tsv}; mv {project_name}.proteinortho-graph {output.graph}; mv {project_name}.info {output.info}; "
        f"mv {project_name}.blast-graph {output.blast_graph}; mv {project_name}.proteinortho-graph.summary {output.graph_summary}; mv {project_name}.proteinortho.html {output.html};"

rule classify_ipr_tfs:
    input:
        ipr_validate=f"{outdir}/.all_iprscan_done",
        proteinortho_tsv=f"{outdir}/proteinortho/{project_name}.proteinortho.tsv"
    output:
        validate=f"{outdir}/.classify_tfs_done",
        tf_proteinortho=expand(f"{outdir}/tf_data/{{tf}}.proteinortho.tsv", tf=tf_families_list)
    conda:
        "workflow/envs/aurtho.yml"
    params:
        ipr_dir=f"{outdir}/interproscan/",
        perc_threshold=config["tf_perc_threshold"],
        tf_families=config["tf_families"]
    benchmark:
        f"{outdir}/benchmarks/classify_ipr_tfs/benchmark.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        f"python workflow/scripts/5_classify_tfs_IPR.py --ipr_dir {params.ipr_dir} --tf_domains resources/tf_ipr_domains.tsv --proteinortho_tsv {input.proteinortho_tsv} --outdir {outdir}/tf_data/ --tf_families {params.tf_families} --perc_threshold {params.perc_threshold}"


# rule transform_proteinortho:
#     input:
#         cds_faa=cds_faa_inputs,  # Use function instead of glob
#         main_tsv=f"{outdir}/proteinortho/{project_name}.proteinortho.tsv",
#         graph=f"{outdir}/proteinortho/{project_name}.proteinortho-graph"
#     output:
#         graph_summary=f"{outdir}/proteinortho/{project_name}.proteinortho-graph.summary",
#         html=f"{outdir}/proteinortho/{project_name}.proteinortho.html",
#         xml=f"{outdir}/proteinortho/{project_name}.proteinortho.xml",
#         validate=f"{outdir}/proteinortho/.snakemake_transform_validate"
#     conda:
#         "workflow/envs/proteinortho.yml"
#     benchmark:
#         f"{outdir}/benchmarks/proteinortho/transform_benchmark.tsv"
#         #repeat("benchmarks/somecommand/{sample}.tsv", 3)
#     shell:
#         "proteinortho_summary.pl '{input.graph}' > '{output.graph_summary}' && "
#         "proteinortho2html.pl {input.main_tsv} {input.cds_faa} > {output.html} && "
#         "proteinortho2xml.pl {input.main_tsv} > {output.xml} && "
#         "touch {output.validate}"

checkpoint grab_proteins:
    input:
        cds_faa=cds_faa_inputs,
        tf_validate=f"{outdir}/.classify_tfs_done",
        tf_tsv=f"{outdir}/tf_data/{{tf}}.proteinortho.tsv"
    output:
        validate=f"{outdir}/tf_data/{{tf}}/.snakemake_grab_proteins_validate"
    params:
        cpus=config["proteinortho_cpu"],
        cog_dir=os.path.join(f"{outdir}/tf_data/{{tf}}", "cogs"),
        min_prot=config["minimum_cog_size"]
    conda:
        "workflow/envs/proteinortho.yml"
    benchmark:
        f"{outdir}/benchmarks/proteinortho/grab_proteins_benchmark_{{tf}}.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        "mkdir -p {params.cog_dir} && "
        "proteinortho_grab_proteins.pl -tofiles={params.cog_dir} "
        "-cpus={params.cpus} -minprot={params.min_prot} "
        "{input.tf_tsv} {input.cds_faa} && "
        "touch {output.validate}"

# def get_orthogroups(wildcards):
#     """Dynamically get list of groups after grab_proteins completes"""
#     return glob_wildcards(f"{outdir}/proteinortho/cog/{project_name}.proteinortho.tsv.OrthoGroup{{group}}.fasta")[1]

rule rename_orthogroup_fastas:
    input:
        f"{outdir}/tf_data/{{tf}}/.snakemake_grab_proteins_validate"
    output:
        validate=touch(f"{outdir}/tf_data/{{tf}}/.cogs.renamed.done")
    shell:
        f"rename -f 's/.proteinortho.tsv./_/g' {outdir}/tf_data/*/cogs/*fasta"

### this has to be tested ! And it was in front of the previous rule before as well. I don't know if it's relevant.
# def get_orthogroups(wildcards):
#     tf_data_dir = os.path.join(outdir, "tf_data")
#     # Use the correct dynamic glob with your timestamped outdir
#     pattern = f"{tf_data_dir}/*/cogs/*.proteinortho.tsv.OrthoGroup*.fasta"
#     orthogroup_paths = glob.glob(pattern)
#     print(f"Orthogroup paths: {orthogroup_paths}")
#     # Extract group wildcard values from filenames
#     groups = [os.path.basename(f).split("OrthoGroup")[1].split(".fasta")[0] for f in orthogroup_paths]
#     #print(f"Groups: {groups}")
#     return groups

rule validate_cog_renaming:
    input:
        f"{outdir}/proteinortho/.snakemake_validate",
        expand(f"{outdir}/tf_data/{{tf}}/.snakemake_grab_proteins_validate", tf=tf_families_list),
        expand(f"{outdir}/tf_data/{{tf}}/.cogs.renamed.done", tf=tf_families_list)
    output:
        f"{outdir}/.all_cogs_renamed"
    shell:
        "touch {output}"

# ---------------------------- #
# ------- Classify TFs ------- #
# ---------------------------- #

rule extract_upstream_regions:
    input:
        f"{outdir}/.gff_db_created",
        tf_validate=f"{outdir}/.classify_tfs_done"
    output:
        validate=touch(f"{outdir}/tf_data/{{tf}}/.ups_done_{{length}}"), # when I use temp(), I need to ensure that the file is created in the shell command with touch
        logs=os.path.join(f"{outdir}/tf_data/{{tf}}/upstream_{{length}}.log")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        stop_cds=config["upstream_stop_cds"],
        length=lambda wildcards: wildcards.length,
        tf_dir=lambda wildcards: f"{outdir}/tf_data/{wildcards.tf}"
    benchmark:
        f"{outdir}/benchmarks/extract_upstream_regions/{{tf}}_{{length}}.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        f"python workflow/scripts/6_fetch_ups.py "
        f"--fna-dir {outdir}/fna/ "
        f"--gff-dir {outdir}/cds_gff/ "
        f"--create-db False "
        f"--tf-dir {params.tf_dir} "
        f"--length {params.length} "
        f"--stop-cds {params.stop_cds} "
        f"--logfile {output.logs}" # && "
        #"touch {output.validate}"

rule validate_extract_ups:
    input:
        expand(f"{outdir}/tf_data/{{tf}}/.ups_done_{{length}}", tf=tf_families_list, length=ups_lengths)
    output:
        f"{outdir}/.extract_upstream_regions_done"
    shell:
        "touch {output}"

# ------------------------ #
# ------- run MEME ------- #
# ------------------------ #
rule run_meme:
    input:
        ups_validate=f"{outdir}/.extract_upstream_regions_done"
        #tf_ups_files=lambda wildcards: f"{outdir}/tf_data/{wildcards.tf}/upstream/{wildcards.meme_target}_upstream_{wildcards.length}.fasta"
    output:
        validate=temp(f"{outdir}/tf_data/{{tf}}/meme/.meme_done_{{motif_length}}_{{meme_mode}}") # when I use temp(), I need to ensure that the file is created in the shell command with touch
    params:
        motif_length=lambda wildcards: wildcards.motif_length,
        meme_mode=lambda wildcards: wildcards.meme_mode,
        tf_dir=lambda wildcards: f"{outdir}/tf_data/{wildcards.tf}"
    benchmark:
        f"{outdir}/benchmarks/meme/{{tf}}_{{motif_length}}_{{meme_mode}}.tsv"
        #repeat("benchmarks/somecommand/{sample}.tsv", 3)
    shell:
        "bash workflow/scripts/meme_loop.sh {params.tf_dir}/upstream/ {params.meme_mode} {params.motif_length} {params.tf_dir} && "
        "touch {output.validate}"

rule validate_meme:
    input:
        expand(f"{outdir}/tf_data/{{tf}}/meme/.meme_done_{{motif_length}}_{{meme_mode}}",
               tf=tf_families_list,
               motif_length=meme_lengths,
               meme_mode=meme_modes)
    output:
        f"{outdir}/.all_meme_done"
    shell:
        "touch {output}"

rule extract_meme_info:
    input:
        f"{outdir}/.all_meme_done"
    output:
        touch(f"{outdir}/.meme_info_extracted")
    conda:
        "workflow/envs/aurtho.yml"
    shell:
        "python workflow/scripts/extract_meme_info.py --tf_data {outdir}/tf_data/ --outpath {outdir}/meme_info/"

# ------------------------------------- #
# ------ Integrate dbCAN and TFs ------ #
# ------------------------------------- #
rule integrate_dbcan_tf:
    input:
        dbcan_done=get_dbcan_output,
        tf_validate=f"{outdir}/.classify_tfs_done"
    output:
        outfile=touch(f"{outdir}/.dbcan_tf_integrated")
    conda:
        "workflow/envs/aurtho.yml"
    params:
        method=config["tf_dbcan_method"],
        extension=config["tf_dbcan_extension"],
        neighbours=config["tf_dbcan_neighbours"]
    benchmark:
        f"{outdir}/benchmarks/integrate_dbcan_tf/benchmark.tsv"
    shell:
        f"python workflow/scripts/integrate_dbcan_tf.py " 
        f"--dbcan_dir {outdir}/dbcan/ "
        f"--gff_db_dir {outdir}/cds_gff/ "
        f"--db_dir resources/dbcan/db/ "
        f"--tf_classification {outdir}/tf_data/TF_classification.tsv "
        f"--method {params.method} "
        f"--extension {params.extension} "
        f"--neighbours {params.neighbours} "
        f"--outdir {outdir}/dbcan_tf/"

# ---------------------------- #
# ------- Benchmarking ------- #
# ---------------------------- #
rule combine_benchmarks:
    input:
        f"{outdir}/.dbcan_tf_integrated"
    output:
        combined=f"{outdir}/benchmarks/combined_benchmarks.tsv"
    conda:
        "workflow/envs/aurtho.yml"
    shell:
        f"python workflow/scripts/combine_benchmarks.py --benchmark_dir {outdir}/benchmarks/ --output_file {output.combined}"

# --------------------------------- #
# ----- Target protein report ----- #
# --------------------------------- #
# target_proteins = get.config("target_proteins", "")
# if len(target_proteins) > 0:
#     rule target_protein_report:
#         input:
#             f"{outdir}/.dbcan_tf_integrated"
#         output:
#             report=f"{outdir}/target_protein_report.tsv"
#         conda:
#             "workflow/envs/aurtho.yml"
#         params:
#             proteins=",".join(target_proteins)
#         shell:
#             f"python workflow/scripts/target_protein_report.py "
#             f"--dbcan_tf_dir {outdir}/dbcan_tf/ "
#             f"--tf_classification {outdir}/tf_data/TF_classification.tsv "
#             f"--target_proteins {params.proteins} "
#             f"--output_file {output.report}"

# ------------------------ #
# ------- Clean up ------- #
# ------------------------ #
# rule cleanup_validation_files:
#     input:
#         f"{outdir}/tf_data/.extract_upstream_regions_done"
#     output:
#         f"{outdir}/tf_data/.cleanup_done"
#     shell:
#         """
#         rm -f results/tf_data/*/.ups_done_*
#         touch {output}
#         """