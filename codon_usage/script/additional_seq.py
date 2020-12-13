additional_seq_dir="input_additional_seq"
unique_count={}
files={
        #'pangolin_blast_orf10':['Pangolin','ORF10'],
        #'bat_blast_orf10':['Bat','ORF10'],
        #'wuhan_blast_orf10':['Wuhan','ORF10'],
        #'wuhan_blast_orfE':['Wuhan','ORFE'],
 	'charite_blast_orfe':['Chairite','ORFE'],
        'charite_blast_orf8':['Chairite','ORF8'],	      
        'charite_blast_orf10':['Chairite','ORF10'],
}

added_seq_fp=additional_seq_dir + "/all_added_seq.fasta"

f_out = open(added_seq_fp, "w")
for file in files:
    in_fp = additional_seq_dir + "/blast_seq/" + file
    print(in_fp)
    current_prefix = files[file][1] + "_" + files[file][0].upper() 
    with open(in_fp, 'r') as f_in:
        for line in f_in:
            line=line.strip()
            line=line.split("\t")
            fasta_title=">" + current_prefix + "_" + line[1]
            fasta_seq=line[-1]
            if fasta_title not in unique_count:
                unique_count[fasta_title]=1
            else:
                print(fasta_title)
                fasta_title = fasta_title+"alt_seq"
            f_out.write(fasta_title + "\n" )
            f_out.write(fasta_seq + "\n" )

f_out.close()
            
            
