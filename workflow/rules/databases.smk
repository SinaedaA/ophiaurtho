import platform
system = platform.system()  # Returns 'Linux', 'Darwin' (macOS), or 'Windows'
print(f"Detected system: {system}")

if system == "Windows":
    raise RuntimeError("This workflow is not supported on Windows. Please use Linux or macOS.")

## database downloads
rule all:
    input:
        "resources/dbcan/db/.snakemake_validate",
        "resources/interproscan/.snakemake_validate"

### dbCAN database download
rule dbcan_db_download:
    output:
        "resources/dbcan/db/.snakemake_validate"
    conda:
        "../envs/dbcan.yml"
    shell:
        "mkdir -p resources/dbcan/db; "
        "run_dbcan database --db_dir resources/dbcan/db && touch {output}" # .snakemake_validate will only be created if the download is successful. If internet fails, rule will fail.

### InterProScan database download
if system == "Darwin":
    print("InterProScan database download on MacOS.")
    rule interproscan_db_download:
        output:
            "resources/interproscan/.snakemake_validate"
        shell:
            "curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz && "
            "curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz.md5 && "
            "md5sum -c interproscan-data-5.75-106.0.tar.gz.md5 && "
            "mkdir -p resources/interproscan/db && "
            "tar -xzf interproscan-data-5.75-106.0.tar.gz -C resources/interproscan/db && "
            "rm interproscan-data-5.75-106.0.tar.gz interproscan-data-5.75-106.0.tar.gz.md5 && "
            "touch {output}" # .snakemake_validate will only be created if the download is successful. If internet fails, rule will fail.
elif system == "Linux":
    print("InterProScan database download skipped on Linux. It will be installed during Interproscan software installation. See InterProScan.md file for instructions.")
    rule install_interproscan:
        output:
            "resources/interproscan/.snakemake_validate"
        shell:
            "mkdir -p resources/interproscan/ && cd resources/;"
            "wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz && "
            "wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz.md5 && "
            "md5sum -c interproscan-5.75-106.0-64-bit.tar.gz.md5 && " # Must return *interproscan-5.75-106.0-64-bit.tar.gz: OK*
            "tar -pxvzf interproscan-5.75-106.0-64-bit.tar.gz && "
            "mv interproscan-5*/* interproscan && "
            "rm -rf interproscan-5.75* && "
            "cd interproscan && "
            "python3 setup.py -f interproscan.properties && cd ../../ && "
            "touch {output}" # .snakemake_validate will only be created if the download is successful. If internet fails, rule will fail
