import platform
system = platform.system()  # Returns 'Linux', 'Darwin' (macOS), or 'Windows'
print(f"Detected system: {system}")

if system == "Windows":
    raise RuntimeError("This workflow is not supported on Windows. Please use Linux or macOS.")

## Determine timestamp and output directory
timestamp = config.get("timestamp")
outdir = f"results/{timestamp}"
print(f"Output directory: {outdir}")

rule make_report:
    input:
        tf_data = directory(f"{outdir}/tf_data/"),
        meme_info = directory(f"{outdir}/meme_info/")
    output:
        report = f"{outdir}/target_protein_report.pdf"
    params:
        target_proteins = config["target_protein"]
    conda:
        "../envs/report.yml"
    shell:
        f"python workflow/scripts/target_protein_report.py --target_proteins {params.target_proteins} --tf_data {input.tf_data} --outdir {outdir} --meme_info {input.meme_info}"