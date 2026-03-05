## Get database paths from environment variables
DBCAN_PATH = os.getenv("DBCAN_DB_PATH", "/data/dbcan")
INTERPROSCAN_PATH = os.getenv("INTERPROSCAN_DB_PATH", "/data/interproscan")
## database downloads
rule all:
    input:
        "/data/dbcan/.snakemake_validate",
        "/data/interproscan/.snakemake_validate",
        f"{INTERPROSCAN_PATH}/.ipr_binary_validate"

### dbCAN database download
rule dbcan_db_download:
    output:
        f"{DBCAN_PATH}/.snakemake_validate"
    conda:
        "../envs/dbcan.yml"
    log:
        f"{DBCAN_PATH}/download.log"
    retries: 3
    shell:
        """
        echo "Downloading dbCAN database to {DBCAN_PATH}..." | tee {log}
        mkdir -p {DBCAN_PATH}
        run_dbcan database --db_dir {DBCAN_PATH}/ --aws_s3 2>&1 | tee -a {log}
        
        # Check no .part files exist
        if ls {DBCAN_PATH}/db/*.part 1> /dev/null 2>&1; then
            echo "ERROR: Incomplete .part files found"
            exit 1
        fi
        touch {output}
        echo "dbCAN database download complete!" | tee -a {log}
        """

### InterProScan database download
rule interproscan_db_download:
    output:
        "/data/interproscan/.snakemake_validate"
    shell:
        """
        # Download InterProScan database
        cd /data/interproscan &&
        curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz && 
        curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz.md5 &&
        md5sum -c interproscan-data-5.75-106.0.tar.gz.md5 &&
        tar -xzf interproscan-data-5.75-106.0.tar.gz &&
        rm interproscan-data-5.75-106.0.tar.gz interproscan-data-5.75-106.0.tar.gz.md5 &&
        touch {output} # .snakemake_validate will only be created if the download is successful. If internet fails, rule will fail.
        """

### Interproscan binary download
rule interproscan_binary_download:
    output:
        f"{INTERPROSCAN_PATH}/.ipr_binary_validate"
    shell:
        """
        mkdir -p {INTERPROSCAN_PATH}
        cd {INTERPROSCAN_PATH}
        curl -O https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz
        curl -O https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz.md5
        md5sum -c interproscan-5.75-106.0-64-bit.tar.gz.md5
        tar -xzf interproscan-5.75-106.0-64-bit.tar.gz
        rm interproscan-5.75-106.0-64-bit.tar.gz*
        touch {output}
        """