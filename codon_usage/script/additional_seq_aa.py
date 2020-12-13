import os

in_fp = "input_additional_seq/all_added_seq.fasta"
out_fp = "input_additional_seq/all_added_seq.faa"

def translate(seq):  
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'', 
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W', 
    } 
    protein ="" 
#    if len(seq)%3 == 0: 
    for i in range(0, len(seq), 3): 
        codon = seq[i:i + 3] 
        if codon in table:
            protein+= table[codon]
        elif len(codon) <3:
            protein+= ""
        else:
            protein+= "X"
            
    return protein 

f_out = open(out_fp, 'w')
with open(in_fp,'r') as f_in:
    current_name = None
    current_seq = ""
    for line in f_in:
        line=line.strip()
        if line.startswith(">"):
            if current_name:
                # break the sequence into multiple lines
                out_aa = translate(current_seq)
                out_aa = [out_aa[i:i+80] for i in range(0,len(out_aa), 80)]
                if len(out_aa)> 0:
                    f_out.write(current_name + "\n")
                    f_out.write("\n".join(out_aa) + "\n")
            current_name= line
            current_seq = ""
            
        else:
            current_seq = current_seq + line
            
            
f_out.close()
