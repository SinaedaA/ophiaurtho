# Use mambaforge as base for faster dependency resolution
FROM condaforge/mambaforge:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

## Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    curl \
    build-essential \
    rename \
    docker.io \
    rsync \
    openjdk-11-jdk \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

## Change SSL certificate config
RUN conda config --set ssl_verify /etc/ssl/certs/ca-certificates.crt

## Check conda SSL configuration and connectivity
RUN echo "=== Checking conda SSL config ===" \
    && conda config --show ssl_verify \
    && echo "=== Checking certificate file exists ===" \
    && ls -la /etc/ssl/certs/ca-certificates.crt \
    && echo "=== Testing conda can download ===" \
    && conda search python --info 2>&1 | head -20

## Set working directory
WORKDIR /app

## Copy environment files first (for better Docker layer caching)
COPY workflow/envs /app/workflow/envs

## Create and activate snakemake environment
RUN mamba create -c conda-forge -c bioconda -n snakemake snakemake -y && \
    mamba clean --all -f -y

## Make snakemake environment available
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

## Copy the rest of the project
COPY . /app

## Ensure snakemake conda environments can be created
RUN mkdir -p /app/.snakemake
RUN conda config --set channel_priority strict

## Set up environment variables
ENV PATH="/opt/conda/envs/snakemake/bin:$PATH"
ENV CONDA_DEFAULT_ENV=snakemake

## Activate snakemake environment by default
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake"]
CMD ["/bin/bash"]