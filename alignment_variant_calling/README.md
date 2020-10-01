# Coronahack2020 alignment_and_variant_calling

This pipeline uses genbank reference file NC_045512.2 for alignment and variant calling of fasta files using minimap2 followed by functional consequence prediction with VEP.

Some hacking required to get the genbank created gff to work with VEP, so some QC may be required.

## Usage:
```
   git clone https://github.com/coronahack2020/alignment_and_variant_calling.git
   
   ./alignment_and_variant_calling/run.sh <fasta_directory> <reference prefix>
```
Reference prefix should correspond to the path to reference FASTA and GFF files minus the .fasta/.gb extension. All files will be output to current working directory.

Final VCF output will be called "combined_variants.vep.annot.bcsq.vcf.gz". This will include functional annotations from VEP and haplotype aware consequences from bcftools csq, plus the following custom INFO fields:

* VAC - the number of samples with each ALT allele (comma separated)
* VariantSamples - the ID of each sample with corresponding ALT alleles (comma separated per ALT, samples are pipe separated)



## FASTA re-headered and split as follows:

    sed s,/,_,g ../phylo/wuhan_bat_pangolin.fasta | \
        perl -wanE 'if (/^\>/){ say $F[0];}else{print;}' \
         > wuhan_bat_pangolin.formatted.fasta
    mkdir split_fasta
    pyfasta split --header "split_fasta/%(seqid)s.fasta" wuhan_bat_pangolin.formatted.fasta

## Variant calling run as follows:

    VAR_REPO=$(realpath  ../../../alignment_and_variant_calling/)
    $VAR_REPO/run.sh split_fasta $VAR_REPO/ref/NC_045512 .fasta