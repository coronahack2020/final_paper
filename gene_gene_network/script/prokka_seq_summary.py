# Work out a summary on the genes identified for each sample/species
import os, re

in_seq_to_subject_fp="sequence_to_subject.txt"
in_clstr_dir="input_individual_clstr"
in_gff_dir = "input_individual_gtf"
in_faa_dir = "input_individual_faa"

out_fp="genes_identified.txt"
out_fp2="genes_count_aa.txt"
sequence_to_subject = {}
genes_identified_per_species= {}

all_files=os.listdir(in_gff_dir)
all_samples = {}
for file in all_files:
    current_file = os.path.join(in_gff_dir, file)
    current_species = file.replace(".gff", "")
    genes_identified_per_species[current_species]={}
    with open(current_file, 'r') as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            else:
                line=line.strip()
                if "ID=" in line:
                    line=line.split()
                    pattern= re.compile(r'ID=(.*?);')
                    match_pattern=re.search(pattern, line[8])
                    if match_pattern:
                        current_gene = match_pattern.group(1)
                        current_subject = line[0]
                        sequence_to_subject[current_gene ]=current_subject
                        genes_identified_per_species[current_species].update({current_gene :[]}) if current_gene not in genes_identified_per_species[current_species] else 0
                        genes_identified_per_species[current_species][current_gene].append(current_subject)
                        all_samples[current_subject]=1
f_out = open(out_fp,'w+')          
for seq in sequence_to_subject:
    f_out.write("%s\t%s\n" %(seq, sequence_to_subject[seq]))    
f_out.close()


# read the aa files and look for sequencces per sample
all_files=os.listdir(in_faa_dir)
sample_aa_dict= {}
for file in all_files:
    current_file = os.path.join(in_faa_dir, file)
    current_species = file.replace(".fsa", "")
    with open(current_file, 'r') as f_in:
        for line in f_in:
            if line.startswith(">") and not(current_species == "WUHANREF"):
                line=line.strip()
                current_prokka_name = line.split(" ")[0].replace(">%s_"%current_species, "")
                current_sample = genes_identified_per_species[current_species][current_prokka_name][0]
                sample_aa_dict.update({current_sample:{}}) if current_sample not in sample_aa_dict else 0
                sample_aa_dict[current_sample][current_prokka_name] = 1

f_out = open(out_fp2, 'w')
f_out.write("Genome\tprokka_count\n")
for sample in sample_aa_dict:
    out_line = "%s\t%s\n" %(sample, len(sample_aa_dict[sample]))
    f_out.write(out_line)
f_out.close()


