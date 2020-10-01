#!/bin/bash 
export VEP=/opt/databiology/apps/ensembl-vep-release-92/vep 
export BCFTOOLS=$(which bcftools) #or provide path if not in PATH
export TABIX=$(which tabix) #or provide path
export SAMTOOLS=$(which samtools) #or provide path
export BOWTIE2=$(which bowtie2) #or provide path
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export REFDIR=$DIR/ref
