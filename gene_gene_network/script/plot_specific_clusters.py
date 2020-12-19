# this script prints the average pident for each pangolin or bat genome vs all SARS-CoV-2
import sys
input_cluster_fp = sys.argv[1]
megablast_result_fp='blast_result/all.fasta.tab'
blastn_result_fp='blast_result/all.fasta.less.stringent.tab'
summarised_results=input_cluster_fp.replace(".csv", "_vs_sars_cov_2_pident.csv")

sequence_check={}
sequence_to_species={}
sequence_to_sample_id={}
blast_similarity={}

# read in sequence ID in the cluster
with open(input_cluster_fp, 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        if line_num>0 and len(line)>0:
            line=line.split(",")
            species=line[4]
            species="WUHAN" if species=="WUHANREF" else species
            sequence_id=line[0]
            sequence_check.update({species:{}}) if species not in sequence_check else 0
            sequence_check[species].update({sequence_id:{}}) if sequence_id not in [species] else 0
            sequence_to_species[sequence_id]=species
            sequence_to_sample_id[sequence_id]=line[3]
# store the percent identy compared to each wuhan or wuhanref genome 
with open(megablast_result_fp, 'r') as f_in:
    for line in f_in:
        line=line.strip()
        line=line.split("\t")
        seq_ids=[line[0], line[1]]
        seq_ids.sort()
        seq_id_a=seq_ids[0]
        seq_id_b=seq_ids[1]
        if (seq_id_a in sequence_to_species) and (seq_id_b in sequence_check['WUHAN']):
            current_species=sequence_to_species[seq_id_a]
            pident=float(line[2])
            sequence_check[current_species][seq_id_a].update({seq_id_b:0}) if seq_id_b not in sequence_check[current_species][seq_id_a] else 0
            # keep the higher pident
            sequence_check[current_species][seq_id_a][seq_id_b] = pident if pident > sequence_check[current_species][seq_id_a][seq_id_b] else sequence_check[current_species][seq_id_a][seq_id_b]

#print(sequence_check)
# print out a csv file with the average percentage identiy for plots
f_out = open(summarised_results, 'w')
f_out.write(",".join(['sample_id', 'species', 'avg_pident', 'num_hits']) + "\n")
for species in sequence_check:
    for sequence_id in sequence_check[species]:
        current_sample=sequence_to_sample_id[sequence_id]
        current_all_pident_list=sequence_check[species][sequence_id].values()
        current_all_pident=sum(current_all_pident_list)/len(current_all_pident_list) if len(current_all_pident_list)> 0 else 0
        f_out.write(",".join([current_sample, species, "%.2f"%current_all_pident, "%d"%len(current_all_pident_list)]) + "\n")
    

f_out.close()
