# this script will use the MCL table to summarise :
# 1) Which genes are found within each sample
# 2) The length of the gene found within each sample
# 3) Print out a graph with x-axis showing the samples, y-axis showing the genes, and fill showing the presence/absence of each gene/length of the gene 

# factors to consider:
# similar genes between sheep/goat/lumpy
# separate plot with just the core genes
# how to deal with genes that are found in both the start/end of the genome

# steps 
# 1) which known genes does each cluster have

import os
working_dir="D:/git/coronahack/PROKKA_blast"
os.chdir(working_dir)
in_mcl_fp = "mcl/nucleotides.csv"
out_summarised_count_fp = "summary_table/gene_count_summary.txt"
if not os.path.exists(os.path.dirname(out_summarised_count_fp)):
    os.mkdir(os.path.dirname(out_summarised_count_fp))
# read in the mcl clustering
organised_data={}
cluster_ref_count = {}
cluster_ref={}
output_genes={}
species_total = {'BAT':215,
                 'PANGOLIN':7,
                 'WUHAN':46,}
# first round go through - label cluster ids and wuhan refs
with open(in_mcl_fp , 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        vals = line.split(",")
        if line_num==0:
            colnames=vals
            nodename_idx =[idx for idx,val in enumerate(colnames) if val == "Node Name"][0]
            sample_idx = [idx for idx,val in enumerate(colnames) if val == "Node sample"][0]
            cluster_idx = [idx for idx,val in enumerate(colnames) if val == "MCL Cluster"][0]
            gene_idx = [idx for idx,val in enumerate(colnames) if val == "Node gene"][0]
        else:
            sample=vals[sample_idx]
            cluster=vals[cluster_idx]
            gene=vals[gene_idx]
            if sample.startswith("WUHANREF"):
                output_genes[cluster]=gene
# second round go through - label prokka gene names with wuhan refs
prokka_to_output_gene={}
with open(in_mcl_fp , 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        vals = line.split(",")
        if line_num==0:
            colnames=vals
            nodename_idx =[idx for idx,val in enumerate(colnames) if val == "Node Name"][0]
            sample_idx = [idx for idx,val in enumerate(colnames) if val == "Node sample"][0]
            cluster_idx = [idx for idx,val in enumerate(colnames) if val == "MCL Cluster"][0]
            gene_idx = [idx for idx,val in enumerate(colnames) if val == "Node gene"][0]
        else:
            sample=vals[sample_idx]
            cluster=vals[cluster_idx]
            gene=vals[gene_idx]
            if cluster in output_genes:
                wuhan_ref_gene = output_genes[cluster]
                if gene != "hypothetical protein":
                    prokka_to_output_gene[gene] = wuhan_ref_gene

# third round through. label cluster with prokka gene name if the cluster/gene name is not labelled with wuhan ref
with open(in_mcl_fp , 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        vals = line.split(",")
        if line_num==0:
            colnames=vals
            nodename_idx =[idx for idx,val in enumerate(colnames) if val == "Node Name"][0]
            sample_idx = [idx for idx,val in enumerate(colnames) if val == "Node sample"][0]
            cluster_idx = [idx for idx,val in enumerate(colnames) if val == "MCL Cluster"][0]
            gene_idx = [idx for idx,val in enumerate(colnames) if val == "Node gene"][0]
        else:
            sample=vals[sample_idx]
            cluster=vals[cluster_idx]
            gene=vals[gene_idx]
            if gene in prokka_to_output_gene and cluster not in output_genes:
                wuhan_ref_gene = prokka_to_output_gene[gene]
                output_genes[cluster] = wuhan_ref_gene 
            if (gene != "hypothetical protein") and (gene not in prokka_to_output_gene):
                prokka_to_output_gene[gene] = gene
                if cluster not in output_genes:
                    output_genes[cluster] = gene 

# forth round going through - work out which wuhan gene does each sample have
output_summary={}
species_count={}
with open(in_mcl_fp , 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        vals = line.split(",")
        if line_num==0:
            colnames=vals
            nodename_idx =[idx for idx,val in enumerate(colnames) if val == "Node Name"][0]
            sample_idx = [idx for idx,val in enumerate(colnames) if val == "Node sample"][0]
            cluster_idx = [idx for idx,val in enumerate(colnames) if val == "MCL Cluster"][0]
            gene_idx = [idx for idx,val in enumerate(colnames) if val == "Node gene"][0]
            species_idx = [idx for idx,val in enumerate(colnames) if val == "Node species"][0]
        else:
            sample=vals[sample_idx]
            cluster=vals[cluster_idx]
            gene=vals[gene_idx]
            species=vals[species_idx]
            species_count.update({species:{}}) if species not in species_count else 0
            # remove extra characters in the sample names
            sample = sample.replace("WUHAN_", "")
            sample = sample.replace("BAT_", "")
            sample = sample.replace("PANGOLIN_", "")
            species_count[species][sample]=1
            if gene in prokka_to_output_gene:
                current_output_gene=prokka_to_output_gene[gene]
            elif cluster in output_genes:
                current_output_gene=output_genes[cluster]
            else:
                continue
            output_summary.update({current_output_gene:{}}) if current_output_gene not in output_summary else 0
            output_summary[current_output_gene].update({species:[]}) if species not in output_summary[current_output_gene] else 0
            output_summary[current_output_gene][species].append(sample) if sample not in output_summary[current_output_gene][species] else 0
            

# Go through the output genes - work out how many samples of each gene
f_out= open(out_summarised_count_fp , 'w')
f_out.write("%s\t%s\t%s\t%s\n"%("gene", "species", "count", "percentage"))
for gene in output_summary:
    for species in output_summary[gene]:
        if species != "WUHANREF":
            current_species_total=species_total[species]
            current_count=len(output_summary[gene][species])
            percentage=current_count/current_species_total * 100
            out_line="%s\t%s\t%d\t%.2f\n"%(gene, species, current_count,percentage)
            f_out.write(out_line)
f_out.close()

