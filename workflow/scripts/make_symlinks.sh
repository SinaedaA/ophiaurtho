#!/bin/bash
assembly_info=$1
genome_dir=$2

cd $genome_dir

IFS=$'\n'
for line in $(cat ../../../$assembly_info); do 
    strain=`echo $line | cut -f1`
    gcf=`echo $line | cut -f2`
    gunzip $gcf/${gcf}_genomic.fna.gz
    cd ../cds_faa/; ln -s ../genomes/$gcf/$strain.cds.faa $strain.cds.faa
    cd ../cds_gff/; ln -s ../genomes/$gcf/$strain.cds.gff $strain.cds.gff
    cd ../fna/; ln -s ../genomes/$gcf/${gcf}_genomic.fna $strain.fna
    cd ../genomes
done