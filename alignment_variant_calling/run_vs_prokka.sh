#!/bin/bash

if [ $# -gt 3 ] || [ $# -lt 2 ]
then
    echo "Usage: $0 fasta_dir ref_prefix <suffix>"
    exit 1
fi

FASTADIR=$1
REF_PREFIX=$2
SUFFIX=$3
if [ "$SUFFIX" == "" ]
then
    SUFFIX=".genome.fasta"
fi
set -euo pipefail
REF_FASTA=${REF_PREFIX}.fasta
REF_GFF=${REF_PREFIX}.gff
if [ ! -f "$REF_FASTA" ]
then
    REF_FASTA=${REF_PREFIX}.fa
    if [ ! -f "$REF_FASTA" ]
    then
        echo "No reference .fasta or .fa found for $REF_PREFIX"
        exit 2
    fi
fi

REF_GFF=${REF_PREFIX}.gff
if [ ! -f "$REF_GFF" ]
then
    echo "$REF_GFF does not exist - exiting"
    exit 3
fi


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source $DIR/config.sh

echo $(date) retrieving minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
export PATH=$PATH:$(realpath minimap2-2.17_x64-linux/)

mkdir -p alignments variants bcfs
$SAMTOOLS faidx $REF_FASTA
echo $(date) aligning .fasta files in $FASTADIR
for FA in $FASTADIR/*.fasta 
do 
    minimap2 --cs -cx asm20 $REF_FASTA $FA > alignments/$(basename $FA $SUFFIX).paf
done

echo $(date) calling variants
for PAF in alignments/*.paf
do 
    sort -k6,6 -k8,8n $PAF  | paftools.js call -l 200 -L 200 -q 30 -s $(basename $PAF .paf) -f $REF_FASTA - > variants/$(basename $PAF .paf).vcf
done

echo $(date) converting to BCF
for VCF in variants/*vcf
do 
    $BCFTOOLS view -O b -o bcfs/$(basename $VCF .vcf).bcf $VCF 
done

echo $(date) indexing BCFs
for BCF in bcfs/*bcf
do 
    $BCFTOOLS index $BCF
done

u=$(expr $(ulimit -n) / 2)
BCFS=(bcfs/*.bcf)
if [ ${#BCFS[@]} -gt $u ]
then
    echo $(date) merging ${#BCFS[@]} BCFs in batches of $u
    mkdir -p merges
    i=0
    for ((n=0; n<${#BCFS[@]}; n+=$u))
    do 
        i=$(expr $i + 1)
        echo $(date) doing merge $i
        SUBSET=${BCFS[@]:$n:$u}
        OUT=merges/merge.${i}.bcf
        bcftools merge -O b -o $OUT $SUBSET
        bcftools index $OUT
    done
    echo $(date) doing final merge
    $BCFTOOLS merge -O z -o combined_variants.vcf.gz merges/*bcf  
else
    echo $(date) merging BCFs
    $BCFTOOLS merge -O z -o combined_variants.vcf.gz bcfs/*bcf  

fi

echo $(date) converting GFF reference for VEP
CSQ_GFF=$(dirname $REF_GFF)/$(basename $REF_GFF .gff).csq.gff
$DIR/convert_prokka_gff.py $REF_GFF > $CSQ_GFF 
grep -v "#" $CSQ_GFF | grep -v '^$' | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip  -c > $CSQ_GFF.gz
$TABIX -p gff $CSQ_GFF.gz
echo $(date) running VEP
#$VEP --fasta $REF_FASTA  --gff $CSQ_GFF.gz -i combined_variants.vcf.gz --force 
$VEP --fasta $REF_FASTA  --gff $CSQ_GFF.gz -i combined_variants.vcf.gz --force --hgvs --compress bgzip --vcf -o combined_variants.vep.vcf.gz

echo $(date) Running bcftools haplotype aware CSQ
$BCFTOOLS csq -O z -f $REF_FASTA -g $CSQ_GFF.gz -p a  combined_variants.vep.vcf.gz > combined_variants.vep.bcsq.vcf.gz
$BCFTOOLS query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ{*}\n]' combined_variants.vep.bcsq.vcf.gz  > combined_variants.bcsq.tsv

echo $(date) Annotating sample call numbers per variant
$DIR/annotate_sample_calls.py combined_variants.vep.bcsq.vcf.gz

echo $(date) getting alignment stats
mkdir -p stats
for PAF in alignments/*paf
do 
    if [ "$(cat $PAF | wc -l)" == "0" ]
    then 
        basename $PAF .paf
    fi 
done > failed_to_align.txt 
for PAF in alignments/*paf
do 
    paftools.js stat $PAF > stats/$(basename $PAF .paf).stats
done

echo $(date) Writing CSV output
$DIR/vcf_to_df.py combined_variants.vep.bcsq.annot.vcf.gz failed_to_align.txt combined_variants.vep.bcsq.annot.csv

echo $(date) Done
