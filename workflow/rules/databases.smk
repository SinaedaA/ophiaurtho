## database downloads
rule all:
    input:
        "resources/dbcan/db/.snakemake_validate",
        "resources/interproscan/db/.snakemake_validate"

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
rule interproscan_db_download:
    output:
        "resources/interproscan/db/.snakemake_validate"
    shell:
        "curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz && "
        "curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz.md5 && "
        "md5sum -c interproscan-data-5.75-106.0.tar.gz.md5 && "
        "mkdir -p resources/interproscan/db && "
        "tar -xzf interproscan-data-5.75-106.0.tar.gz -C resources/interproscan/db && "
        "rm interproscan-data-5.75-106.0.tar.gz interproscan-data-5.75-106.0.tar.gz.md5 && "
        "touch {output}" # .snakemake_validate will only be created if the download is successful. If internet fails, rule will fail.