# organise a sequence to sample conversion
import os, re

in_dir = "input_individual_gtf"
out_fp="sequence_to_subject.txt"
sequence_to_subject = {}

all_files=os.listdir(in_dir)

for file in all_files:
    current_file = os.path.join(in_dir, file)
    current_species = file.replace(".gff", "")
    with open(current_file, 'r') as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            else:
                line=line.strip()
                if "ID=" in line:
                    vals=line.split("\t")
                    pattern_ID= re.compile(r'ID=(.*?);')
                    pattern_gene= re.compile(r'.*?;product=(.*?)$')
                    match_pattern_ID=re.search(pattern_ID, vals[-1])
                    match_pattern_gene=re.search(pattern_gene, vals[-1])
                    start=int(vals[3])
                    end=int(vals[4])
                    if match_pattern_ID and match_pattern_gene:
                        sequence_to_subject[match_pattern_ID.group(1)]=[vals[0], current_species, match_pattern_gene.group(1),str(end-start)]
                        
f_out = open(out_fp,'w+')          
f_out.write("seq_id\tsample\tspecies\tgene\tgeneLength\n")
for seq in sequence_to_subject:
    f_out.write("%s\t%s\n" %(seq, "\t".join(sequence_to_subject[seq])))
    
    
f_out.close()
