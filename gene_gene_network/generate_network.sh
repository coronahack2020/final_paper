#### REQUIREMENT ####
# This script requires blast+ (version 2.7.1) and python (version 3.4.3)

# this script generate the network analysis.gml files for Graphia
PROKKA_DIR=../host-data

INDIVIDUAL_FASTA_DIR=input_individual_fasta
INDIVIDUAL_FAA_DIR=input_individual_faa
INDIVIDUAL_TSV_DIR=input_individual_tsv
INDIVIDUAL_GTF_DIR=input_individual_gtf
ADDITIONAL_SEQ=input_additional_seq
COLLATED_SEQ=input_collated
BLAST_RESULT_DIR=blast_result
REF_GENES=input_ref_genes # the files in this folder are curated from Ensembl v100 ASM985889v3 
 
#module add roslin/blast+/2.7.1
#module add python/3.4.3


mkdir -p $INDIVIDUAL_FASTA_DIR $INDIVIDUAL_FAA_DIR $INDIVIDUAL_TSV_DIR $INDIVIDUAL_GTF_DIR $ADDITIONAL_SEQ $ADDITIONAL_SEQ/blast_seq $COLLATED_SEQ $BLAST_RESULT_DIR

#### Organise genes identified by Prokka
# get the sample annotation lines from the ffn files
cp $PROKKA_DIR/bat/bat_prokka/*.ffn $INDIVIDUAL_FASTA_DIR/BAT.fasta
cp $PROKKA_DIR/pangolin/pangolin_prokka/*_prokka.ffn $INDIVIDUAL_FASTA_DIR/PANGOLIN.fasta
cp $PROKKA_DIR/wuhan/wuhan_prokka/*_prokka.ffn $INDIVIDUAL_FASTA_DIR/WUHAN.fasta

# get the sample annotation lines from the faa files
cp $PROKKA_DIR/bat/bat_prokka/*.faa $INDIVIDUAL_FAA_DIR/BAT.fsa
cp $PROKKA_DIR/pangolin/pangolin_prokka/*_prokka.faa $INDIVIDUAL_FAA_DIR/PANGOLIN.fsa
cp $PROKKA_DIR/wuhan/wuhan_prokka/*_prokka.faa $INDIVIDUAL_FAA_DIR/WUHAN.fsa

# get the annotation lines from the fnn files
cp $PROKKA_DIR/bat/bat_prokka/*.tsv $INDIVIDUAL_TSV_DIR/BAT.tsv
cp $PROKKA_DIR/pangolin/pangolin_prokka/*_prokka.tsv $INDIVIDUAL_TSV_DIR/PANGOLIN.tsv
cp $PROKKA_DIR/wuhan/wuhan_prokka/*_prokka.tsv $INDIVIDUAL_TSV_DIR/WUHAN.tsv

# get the sample annotation lines from the gff files
cp $PROKKA_DIR/bat/bat_prokka/*.gff $INDIVIDUAL_GTF_DIR/BAT.gff
cp $PROKKA_DIR/pangolin/pangolin_prokka/*_prokka.gff $INDIVIDUAL_GTF_DIR/PANGOLIN.gff
cp $PROKKA_DIR/wuhan/wuhan_prokka/*_prokka.gff $INDIVIDUAL_GTF_DIR/WUHAN.gff

#### Add the Ensembl v100 ASM985889v3 to the folders
cp $REF_GENES/WUHANREF.fasta $INDIVIDUAL_FASTA_DIR
cp $REF_GENES/WUHANREF.fsa $INDIVIDUAL_FAA_DIR

# change the annotation to include the host species
for file in $INDIVIDUAL_FASTA_DIR/*.fasta; do
replacement_string=$(basename $file)
replacement_string=${replacement_string//.fasta/}
replacement_string=${replacement_string//_/}
replacement_string=">${replacement_string}_"
sed -i "s/>/$replacement_string/g" $file
done

# change the annotation to include the host species
for file in $INDIVIDUAL_FAA_DIR/*.fsa; do
replacement_string=$(basename $file)
replacement_string=${replacement_string//.fsa/}
replacement_string=${replacement_string//_/}
replacement_string=">${replacement_string}_"
sed -i "s/>/$replacement_string/g" $file
done

#### Organise genes identified by Prokka
# Add additional sequences identified through blast 
# (adding orf10 and orfE - orf8 were not added because the three identified were already identified by prokka) to the data
# Convert the blast output into fasta
cp $PROKKA_DIR/*/*blast* $ADDITIONAL_SEQ/blast_seq
python script/additional_seq.py
python script/additional_seq_aa.py

#### merge the sequences together into one file
cat $INDIVIDUAL_FASTA_DIR/*.fasta $ADDITIONAL_SEQ/all_added_seq.fasta > $COLLATED_SEQ/all.fasta
cat $INDIVIDUAL_FAA_DIR/*.fsa $ADDITIONAL_SEQ/all_added_seq.faa > $COLLATED_SEQ/all_aa.fsa


#### Run gene-gene BLAST analysis 
# using the nucleotide data
MERGED_FSA=$COLLATED_SEQ/all.fasta
BLAST_OUTPUT=blast_result/all.fasta.tab
rm ${MERGED_FSA}.* # remove existing blastdb
makeblastdb -in $MERGED_FSA -parse_seqids -dbtype nucl

blastn -db "${MERGED_FSA}" \
-query $MERGED_FSA \
-out $BLAST_OUTPUT \
-outfmt '6 qseqid sseqid pident sstrand evalue qlen slen length qcovs mismatch gapopen bitscore'

# using the amino acid data
MERGED_FSA=$COLLATED_SEQ/all_aa.fsa
BLAST_OUTPUT=blast_result/all_aa.faa.tab
rm ${MERGED_FSA}.* # remove existing blastdb
makeblastdb -in $MERGED_FSA -parse_seqids -dbtype prot

blastp -db "${MERGED_FSA}" \
-query $MERGED_FSA \
-out $BLAST_OUTPUT \
-outfmt '6 qseqid sseqid pident sstrand evalue qlen slen length qcovs mismatch gapopen bitscore'


#### Make network graph
python script/sequence_to_sample.py
python script/prokka_seq_summary.py
python script/blast_to_self_summarise_nuc.py
python script/blast_to_self_summarise_aa.py

#### Open the all_nuc.gml and all_aa.gml for nucleotide and amino acid network graphs respectively
# This is done using Graphia (https://graphia.app/)
# Add the following transforms:
# Remove edges
# < 60 blastscore
# < 80 qcovs
# < 80 pident
# Remove components
# < 5
# MCL cluster with 2.0 Granularity
















