import platform
system = platform.system()  # Returns 'Linux', 'Darwin' (macOS), or 'Windows'
print(f"Detected system: {system}")

if system == "Windows":
    raise RuntimeError("This workflow is not supported on Windows. Please use Linux or macOS.")

## Determine timestamp and output directory
timestamp = config.get("timestamp")
outdir = f"results/{timestamp}"
print(f"Output directory: {outdir}")

rule grab_all_proteins:
    output:
        validate=f"{outdir}/.created_all_cogs"
    params:
        min_prot=config["minimum_cog_size"],
        cpus=config["proteinortho_cpu"],
        project_name=config["proteinortho_project"]
    conda:
        "../envs/proteinortho.yml"
    shell:
        f"mkdir -p {outdir}/proteinortho/cogs && "
        f"proteinortho_grab_proteins.pl -tofiles={outdir}/proteinortho/cogs/ "
        f"-cpus={params.cpus} -minprot={params.min_prot} "
        f"{outdir}/proteinortho/{params.project_name}.proteinortho.tsv {outdir}/cds_faa/*.cds.faa && "
        f"touch {output.validate} && "
        f"rename 's/\.proteinortho.tsv\./_/g' {outdir}/proteinortho/cogs/*fasta"
