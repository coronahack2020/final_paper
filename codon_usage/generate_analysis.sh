# this script copy over files from the gene-gene network analysis. 
# The scripts in the gene-gene network analysis will need to be generated first.
ORGANISED_ALL_FNN="../gene_gene_network/input_collated/all.fasta"
ORGANISED_ALL_FAA="../gene_gene_network/input_collated/all_aa.fsa"
CHAIRITE_PROKKA_FNN_ORIGIN="../host-data/charite/prokka_out/human_covid_charite_prokka.ffn"
CHAIRITE_PROKKA_FAA_ORIGIN="../host-data/charite/prokka_out/human_covid_charite_prokka.faa"
ORGANISED_SAMPLE_TO_SUBJECT_ORIGIN="../gene_gene_network/sequence_to_subject.txt"
BLAST_SEQ_ADD_ORIGIN=""


# set the fp for input/output files in this directory
MODIFIED_REF_CDNA="input/blast/modified_ref_cdna.fasta"
ALL_SEQ_FNN="input/blast/all.fasta"
ALL_SEQ_FAA="input/blast/all_aa.fsa"
BLAST_OUTPUT="input/blast/blast_result.tab"
CHAIRITE_PROKKA_FNN="input/blast/additional_prokka/human_covid_charite_prokka.ffn"
CHAIRITE_PROKKA_FAA="input/blast/additional_prokka/human_covid_charite_prokka.faa"
ORGANISED_SAMPLE_TO_SUBJECT="input/blast/sequence_to_subject.txt"
ADDITIONAL_SEQ_DIR='input_additional_seq'
ADDITIONAL_SEQ_FNN="input_additional_seq/all_added_seq.fasta"
ADDITIONAL_SEQ_FAA="input_additional_seq/all_added_seq.faa"


mkdir -p "input" "input/blast" "input/blast/additional_prokka" $ADDITIONAL_SEQ_DIR
cp $ORGANISED_ALL_FNN $ALL_SEQ_FNN
cp $ORGANISED_ALL_FAA $ALL_SEQ_FAA
cp $CHAIRITE_PROKKA_FNN_ORIGIN $CHAIRITE_PROKKA_FNN
cp $CHAIRITE_PROKKA_FAA_ORIGIN $CHAIRITE_PROKKA_FAA
sed -i "s/>/>CHAIRITE_/g" $CHAIRITE_PROKKA_FAA
sed -i "s/>/>CHAIRITE_/g" $CHAIRITE_PROKKA_FNN
cp -r "../gene_gene_network/input_additional_seq/blast_seq" $ADDITIONAL_SEQ_DIR
cp modified_ref_cdna.fasta input/blast
# add the additional blast sequencce
python script/additional_seq.py
python script/additional_seq_aa.py
cat $ALL_SEQ_FNN $CHAIRITE_PROKKA_FNN $ADDITIONAL_SEQ_FNN > tmp
mv tmp $ALL_SEQ_FNN
cat $ALL_SEQ_FAA $CHAIRITE_PROKKA_FAA $ADDITIONAL_SEQ_FAA > tmp
mv tmp $ALL_SEQ_FAA
cp $ORGANISED_SAMPLE_TO_SUBJECT_ORIGIN $ORGANISED_SAMPLE_TO_SUBJECT

module unload python
module add roslin/blast+/2.7.1
module add roslin/bedtools/2.29.2
module add python/3.4.3


#### Run (gene)-to-(modified reference cDNA) BLAST analysis 
# using the nucleotide data
rm ${MODIFIED_REF_CDNA}.* # remove existing blastdb
rm $BLAST_OUTPUT
makeblastdb -in $MODIFIED_REF_CDNA -parse_seqids -dbtype nucl
makeblastdb -in $ALL_SEQ_FNN -parse_seqids -dbtype nucl

blastn -db "${MODIFIED_REF_CDNA}" \
	-query $ALL_SEQ_FNN \
	-out $BLAST_OUTPUT \
	-outfmt '6 qseqid sseqid pident sstrand evalue qlen slen length qcovs mismatch gapopen bitscore' 
#### Select for PROKKA ID that are also in the aa file (i.e. retaining only the cDNA that gets translated)
python script/find_best_blast_hit.py
Rscript script/plot_rscu_pca.r # please make sure you have the required libraries installed

