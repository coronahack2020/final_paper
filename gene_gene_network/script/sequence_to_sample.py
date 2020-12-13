import os, re

in_dir = "input_individual_gtf"
out_fp="sequence_to_subject.txt"
chairite_conversion="addition_sample_info/chairite_conversion.txt"
sequence_to_subject = {}
chairite_id_conversion=[]
all_files=os.listdir(in_dir)

for file in all_files:
    current_file = os.path.join(in_dir, file)
    current_species = file.replace(".gff", "")
    if current_species == "CHAIRITE":
        count_idx=0
        with open(chairite_conversion, 'r') as f_in:
            for line in f_in:
                line=line.strip()
                chairite_id_conversion.append(line)

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
                    sample_name = vals[0]
                    if current_species == "CHAIRITE":
                        chairite_idx = int(sample_name.split("_")[-1])-1
                        sample_name = chairite_id_conversion[chairite_idx]
                    if match_pattern_ID and match_pattern_gene:
                        sequence_to_subject[match_pattern_ID.group(1)]=[sample_name, current_species, match_pattern_gene.group(1),str(end-start)]

f_out = open(out_fp,'w+')
f_out.write("seq_id\tsample\tspecies\tgene\tgeneLength\n")
for seq in sequence_to_subject:
    f_out.write("%s\t%s\n" %(seq, "\t".join(sequence_to_subject[seq])))


f_out.close()

