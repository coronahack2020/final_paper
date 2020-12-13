import os, numpy
in_blast_result = "input/blast/blast_result.tab"
in_prokka_aa = "input/blast/all_aa.fsa"
in_prokka_nuc = "input/blast/all.fasta"
in_sequence_to_subject = "input/blast/sequence_to_subject.txt"
out_keep_sequence = "input/blast/keep_prokka_annotation.csv"
out_fasta_by_host_dir = "input/fasta_by_host"
out_fasta_by_gene_dir = "input/fasta_by_gene"
out_orf1ab_merge = "output/summary/orf1ab_length.csv"
out_sample_gene_count = "output/summary/sample_gene_count.csv"
out_sample_gene_length = "output/summary/sample_gene_length.csv"
out_codon_usage = "input/codon_usage_best_gene.csv"

fasta_line_len = 80

# these are incomplete sequeces found by blast
bad_seq_ids = ["ORFE_WUHAN_hCoV-19/Wuhan/WH05/2020alt_seq",
               "ORFE_WUHAN_hCoV-19/Wuhan/WH05/2020", 
               'ORFE_CHAIRITE_BetaCoV/Rhein-Kreis-Neuss/MD572006610/2020',
               "ORF10_CHAIRITE_BetaCoV/Germany/ChVir3214/2020",
               "ORF10_BetaCoV/Germany/ChVir3193/2020", 
               "ORF10_CHAIRITE_BetaCoV/Germany/ChVir3214/2020alt_seq", 
        ]
  
prokka_aa_presence={}
if not os.path.exists(out_fasta_by_host_dir):
    os.mkdir(out_fasta_by_host_dir)
if not os.path.exists(out_fasta_by_gene_dir):
    os.mkdir(out_fasta_by_gene_dir)
if not os.path.exists("output"):
    os.mkdir("output")
if not os.path.exists(os.path.dirname(out_orf1ab_merge)):
    os.mkdir(os.path.dirname(out_orf1ab_merge))

# Functions
def codon_count(f_seq):
    f_codon_summary={}
    f_seq = f_seq.upper()
    f_seq=[f_seq[i:i+3] for i in range(0,len(f_seq), 3)]
    for f_codon in f_seq:
        if len(f_codon) == 3 and 'N' not in f_codon and '-' not in f_codon:
            f_codon_summary.update({f_codon:0}) if f_codon not in f_codon_summary else 0
            f_codon_summary[f_codon]+=1
    return(f_codon_summary)


# read in the numbers in prokka aa to make sure that it is a gene that can be transcribed
with open(in_prokka_aa, 'r') as f_in:
    for line in f_in:
        if line.startswith(">"):
            line=line.strip()
            line=line.split(" ")
            prokka_id=line[0].replace(">", "")
            prokka_aa_presence[prokka_id]=1
# blast output header 
# qseqid sseqid pident sstrand evalue qlen slen length qcovs mismatch gapopen bitscore

# Read in sequence to subject table
prokka_to_sample={}
sample_to_species={'BetaCoV/Germany/ChVir3193/2020alt_seq':"CHAIRITE"}
species_to_sample={}
sarscov2_genes1 = ['ORF1ab','E', 'M', 'N', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF10', 'S', 'ORF1a']
sarscov2_genes2 = ['ORF1a', 'ORF1b'] # these need to be merged together to form ORF1ab

with open(in_sequence_to_subject, 'r') as f_in:
    for line in f_in:
        line=line.strip()
        line=line.split("\t")
        prokka_id = line[0]
        sample_id = line[1]
        species = line[2]
        prokka_to_sample[prokka_id] = sample_id
        sample_to_species[sample_id] = species
        species_to_sample.update({species:{}}) if species not in species_to_sample else 0
        species_to_sample[species].update({sample_id:1}) if sample_id not in species_to_sample[species] else 0
        
# Go through blast result to find the longest sequencce that matches to the reference cdna
blast_output_keep={}
fasta_id_to_species = {}
fasta_id_to_gene = {}
with open(in_blast_result, 'r') as f_in:
    for line in f_in:
        line = line.strip()
        line = line.split("\t")
        full_id = line[0]
        sample_id_split = line[0].split("_")
        current_length=int(line[7])
        current_gene=line[1]
    
        sample_id_0 = sample_id_split[0]
        sample_id_1 = "_".join(sample_id_split[1:len(sample_id_split)])
        if sample_id_0 == "WUHANREF" or full_id in bad_seq_ids: # ignore the sample if it is wuhan ref or one of the bad blast sequences
            continue
        # if this is a prokka sequence (rather than blast), check if it is in the prokka aa file for protein coding
        if sample_id_0 in ['BAT', 'WUHAN', 'PANGOLIN', 'CHAIRITE'] :
            if full_id in prokka_aa_presence:
                current_sample= prokka_to_sample[sample_id_1]
                blast_output_keep.update({current_sample:{}}) if current_sample not in blast_output_keep else 0
                if current_gene not in blast_output_keep[current_sample]:
                    blast_output_keep[current_sample][current_gene]=[current_length, full_id, sample_id_0, sample_id_1] 
                else:
                    if current_length > blast_output_keep[current_sample][current_gene][0]:
                        blast_output_keep[current_sample][current_gene]=[current_length, full_id, sample_id_0, sample_id_1] 
            fasta_id_to_species[full_id]=sample_id_0
        # blast defined genes
        elif not full_id in bad_seq_ids  and "alt_seq" not in full_id:
            sample_id_1_split = sample_id_1.split("_")
            sample_id = "_".join(sample_id_1_split[1:len(sample_id_1_split)])
            dataset_name = sample_id_1_split[0]
            current_sample = sample_id
            blast_output_keep.update({current_sample:{}}) if current_sample not in blast_output_keep else 0
            blast_output_keep[current_sample][current_gene]=[current_length, full_id, dataset_name, sample_id] 
            fasta_id_to_species[full_id]=dataset_name

        fasta_id_to_gene[full_id] = current_gene
# By this point, the only fasta_id in here would be the longest matching sequence for each gene
# Print out an organised table showing the sequence - record the fasta_id for ORF1a annd ORF1b so they can be merged as ORF1ab
keeping_fastaid={}
best_fastaid_per_gene_per_sample={}
f_out = open(out_keep_sequence, 'w')
f_out.write("fasta_id,dataset_id,sample_id,sample_name,gene\n")
for current_sample in blast_output_keep:
    for current_gene in blast_output_keep[current_sample]:
        fasta_id=blast_output_keep[current_sample][current_gene][1]
        dataset_id=blast_output_keep[current_sample][current_gene][2]
        sample_id=blast_output_keep[current_sample][current_gene][3]
        # record the best fastaid corresponding to each gene
        best_fastaid_per_gene_per_sample.update({current_sample:{}}) if current_sample not in best_fastaid_per_gene_per_sample else 0
        best_fastaid_per_gene_per_sample[current_sample].update({current_gene:fasta_id}) 
        keeping_fastaid[fasta_id]=1
        f_out.write(','.join([fasta_id,dataset_id,sample_id,current_sample,current_gene]) + "\n")

f_out.close()
# Go through the fasta file and print out only the sequences where the fasta_id in in the file
# Read in sequences that are preset in keeping_fastaid
fasta_id_to_seq={}
with open(in_prokka_nuc, 'r') as f_in:
    write_lines = False
    for line in f_in:
        line=line.strip()
        if line.startswith(">"):
            vals=line.split(" ")
            current_fasta_id = vals[0].replace(">","")
            if current_fasta_id in keeping_fastaid:
                current_species = fasta_id_to_species[current_fasta_id]
                fasta_id_to_seq[current_fasta_id]=''
                write_lines = True
            else:
                write_lines = False
                current_species = None
        elif write_lines:
            fasta_id_to_seq[current_fasta_id] = fasta_id_to_seq[current_fasta_id] + line
        else:
            continue

# Go through each species, then each sample in each species
# print out in fasta 
# 0- new fasta name that comprised of gene name followed by sample_name
# 1- the fasta sequence in the relavent output
# 2- total number of gene in that sample
# also, merge ORF1a and ORF1b and record the sequences and codon count for each gene along the way
# print out in ORF1ab_conversion whether if both ORF1a and ORF1b exist
gene_len = {}
codon_summary={}
new_fastaid_by_host={}
new_fastaid_by_gene={}
new_fastaid_seq={}
current_gene_count={}
f_out_sample_gene_count = open(out_sample_gene_count, 'w')
f_out_sample_gene_count.write("%s,%s\n"%("sample_name", "gene_count"))

# check through the quality of the genes, filter it and determine if it would be kept for the codon usage dataset
# the following block is for recording gene sequence, gene length 
for current_sample in best_fastaid_per_gene_per_sample:
    if current_sample in sample_to_species:
        current_species=sample_to_species[current_sample]
    else:
        # the charite sample names are sometimes too long so their sample name conversion is sometimes incorrect. Any samples that do not have matching sample names are charite
        current_species="CHAIRITE"
    # print out the rest of the genes 
    for current_gene in sarscov2_genes1:
        # if the gene is ORF1ab, check that both ORF1a annd ORF1b exists for this sample
        current_gene_check=[current_gene] if "ORF1ab" != current_gene else   ['ORF1a','ORF1b']
        best_fastaid_per_gene_per_sample_checked = [val in best_fastaid_per_gene_per_sample[current_sample] for val in current_gene_check]
        best_fastaid_per_gene_per_sample_checked = sum(best_fastaid_per_gene_per_sample_checked)
        if best_fastaid_per_gene_per_sample_checked == len(current_gene_check):
            # find the corresponding fasta_ids and check if they're both in "fasta_id_to_seq
            current_fasta_id = [best_fastaid_per_gene_per_sample[current_sample][val] for val in current_gene_check]
            current_fasta_id_checked = [val in fasta_id_to_seq for val in current_fasta_id]
            current_fasta_id_checked = sum(current_fasta_id_checked)
            if  current_fasta_id_checked == len(current_gene_check):
                new_fasta_id="%s_%s"%(current_gene,current_sample)
                sequence=[fasta_id_to_seq[val] for val in  current_fasta_id]
                # merge the sequence if it is ORF1ab
                sequence=sequence[0] if len(sequence) ==1 else sequence[0] + sequence[1]
                codon_summary[new_fasta_id]=codon_count(sequence) # record codon usage
                gene_len.update({current_gene:{}}) if current_gene not in gene_len else 0 # record gene length
                gene_len[current_gene][new_fasta_id]=len(sequence)   # record gene length
                current_gene_count[current_gene]=1
                # add the new fastaid to relavent dictionsaies
                new_fastaid_by_host.update({current_species:[]}) if current_species not in new_fastaid_by_host else 0
                new_fastaid_by_host[current_species].append(new_fasta_id)
                new_fastaid_by_gene.update({current_gene:[]}) if current_gene not in new_fastaid_by_gene else 0
                new_fastaid_by_gene[current_gene].append(new_fasta_id)
                new_fastaid_seq[new_fasta_id]=sequence
    # print out a summary on how many genes are found in this sample
    f_out_sample_gene_count.write("%s,%d\n"%(current_sample, len(current_gene_count)))
f_out_sample_gene_count.close()

# Go through the gene lengths and print out genes with both passed and failed length (and label the ones that failed)
# if the standard deviation of the sequence length is > 1000, use the criteria keep_length > mean - 0.5*std
# otherwise, use the criteria keep_length > mean - 2*std
new_fastaid_length_passed={}
f_out = open(out_sample_gene_length, 'w')
f_out.write(",".join(['gene_sample_name','gene','genome', 'length', 'qc_length_check']) + "\n")
for gene in gene_len:
    # find the mean, std and the threshold for the length qc
    current_values = list(gene_len[gene].values())    
    current_mean = numpy.mean(current_values)
    current_std = numpy.std(current_values)
    length_qc_threshold = current_mean - 2*current_std if current_std<1000 else current_mean - 0.5*current_std
    length_qc_threshold = 0 if current_std ==0 else length_qc_threshold # no filtering if the standard deviation is less than 2
    for fasta_id in gene_len[gene]:
        length_qc = "passed" if  gene_len[gene][fasta_id] >= length_qc_threshold else "failed"
        # add to the dictionary recording those that passed the length qc
        if length_qc == "passed":
            new_fastaid_length_passed[fasta_id] = 1
        fasta_id_split=fasta_id.split("_")
        sample_id="_".join(fasta_id_split[1:len(fasta_id_split)])
        out_line=[fasta_id, gene, sample_id, str(gene_len[gene][fasta_id]), length_qc]
        f_out.write(",".join(out_line) + "\n")
f_out.close()


# print out fasta by species
for current_species in new_fastaid_by_host:
    current_species_short ='human' if current_species in ['CHAIRITE', 'WUHAN'] else current_species.lower()
    current_out_fp = out_fasta_by_host_dir + "/" + "%s_virus_filtered.fasta" %current_species_short
    f_out = open(current_out_fp, 'w')
    for current_fastaid in new_fastaid_by_host[current_species]:
        if current_fastaid in new_fastaid_length_passed:
            current_seq = new_fastaid_seq[current_fastaid]
            out_current_fastaid = ">" + current_fastaid
            out_current_seq=[current_seq[i:i+fasta_line_len] for i in range(0,len(current_seq), fasta_line_len)]
            f_out.write(out_current_fastaid + "\n")
            f_out.write("%s\n"%("\n".join(out_current_seq)))
    f_out.close()

# print out fasta by gene
for current_gene in new_fastaid_by_gene:
    f_out = open(os.path.join(out_fasta_by_gene_dir, current_gene + ".fasta"), 'w')
    for current_fastaid in new_fastaid_by_gene[current_gene]:
        if current_fastaid in new_fastaid_length_passed:
            current_seq = new_fastaid_seq[current_fastaid]
            out_current_fastaid = ">" + current_fastaid
            out_current_seq=[current_seq[i:i+fasta_line_len] for i in range(0,len(current_seq), fasta_line_len)]
            f_out.write(out_current_fastaid + "\n")
            f_out.write("%s\n"%("\n".join(out_current_seq)))      
    f_out.close()

# print out codon usage summary
all_codon = {}
# get a list of all codon so the print out can be ordered
for current_fastaid in codon_summary:
    for codon in codon_summary[current_fastaid]:
        all_codon[codon] =1
all_codon = list(all_codon.keys())
all_codon = sorted(all_codon)
f_out = open(out_codon_usage, 'w')
f_out.write("%s,%s\n"%("gene_sample_name",",".join(all_codon)))
for current_fastaid in codon_summary:
    if current_fastaid in new_fastaid_length_passed:
        out_line=[current_fastaid]
        for codon in all_codon:
            if codon in codon_summary[current_fastaid]:
                out_line.append(codon_summary[current_fastaid][codon])
            else:
                out_line.append(0)
        f_out.write("%s\n"%(",".join([str(val) for val in out_line])))
f_out.close()


