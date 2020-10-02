import os


in_fp = "blast_result/all.fasta.tab"
in_seq_to_subject_fp="sequence_to_subject.txt"

gene_annotation_dir="input_individual_tsv"

file_suffix="_nuc"
# out files
out_flattened_df_edges = "all_edges_best_blast_score_per_sample%s.csv" %file_suffix
out_flattened_df_nodes="all_nodes_best_blast_score_per_sample%s.csv" %file_suffix
out_gml = "all%s.gml" %file_suffix


    
# print out graphia network
all_nodes_id2name={}
all_nodes_name2id = {}
all_edges ={}
best_blast_score = {}
best_blast_score_node = {}
sample_species={}
seq_details={}
allsub={}
with open(in_seq_to_subject_fp, 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        line=line.split("\t")
        if line_num ==0:
            header = line
        else:
            for header_idx, current_header in enumerate(header):
                seq_details.update({line[0]:{}}) if line[0] not in seq_details else 0
                seq_details[line[0]][current_header]= line[header_idx].replace("/","_")
                allsub[line[header_idx].replace("/","_")] = 1
genomename_to_species={
        'DBBAT': 'DBBAT',
        'VIPRBat': 'VIPRBat',
        'HumanCOVIDCovGLue': 'HumanCOVIDCovGLue',
        'HumanCOVIDCharite': 'HumanCOVIDCharite',
        'DBHumanCOVID': 'DBHumanCOVID',
        'DBHumanMERS': 'DBHumanMERS',
        'NCBIPangolin': 'NCBIPangolin',
        'DBPangolin': 'DBPangolin',
        'BAT': "BAT",
        'PANGOLIN': 'PANGOLIN',
        'WUHAN': 'WUHAN',

        }
genomename_to_species={
        'DBBAT': 'BAT',
        'VIPRBat': 'BAT',
        'HumanCOVIDCovGLue': 'HUMAN_COVID',
        'HumanCOVIDCharite': 'HUMAN_COVID',
        'DBHumanCOVID': 'HUMAN_COVID',
        'DBHumanMERS': 'HUMAN_MERS',
        'NCBIPangolin': 'PANGOLIN',
        'DBPangolin': 'PANGOLIN',
        'BAT': "BAT",
        'PANGOLIN': 'PANGOLIN',
        'WUHAN': 'WUHAN',
        'WUHANREF': 'WUHANREF',

        }
sample_from_each_group={
        'DBBAT':[],
        'VIPRBat': [],
        'HumanCOVIDCovGLue': [],
        'HumanCOVIDCharite': [],  
        'DBHumanCOVID': [],  
        'DBHumanMERS': [],  
        'NCBIPangolin': [],  
        'DBPangolin': [],   
        'BAT': [],  
        'PANGOLIN': [],  
        'WUHAN': [],  
        }
gene_name_to_wuhanref_gene_name={
        'ORFE':'E',
        'Protein 7a': 'ORF7a',
        'Spike glycoprotein': 'S',
        'Membrane protein':'M',
        'Protein 3a':'ORF3a',
        'Nucleoprotein':'N',
        'Enveolope small membrane protein':'E',        
        }
sample_to_species = {}
all_species=[]
# read in gene annotation
all_gene_annotation={}
all_gene_annotation_fp = os.listdir(gene_annotation_dir)
for file in all_gene_annotation_fp:
    current_species = file.replace(".tsv","").replace("_","")
    all_gene_annotation[current_species]={}
    with open(os.path.join(gene_annotation_dir, file), 'r') as f_in:
        for line_num, line in enumerate(f_in):
            if line_num > 0:
                line = line.strip()
                line = line.split("\t")
                all_gene_annotation[current_species][line[0]] = line[-1]
all_gene_annotation.keys()   

filtered_blast_results={}
node_id=1
node_name2id={}
node_id2name={}
node_annotation={}

#   function for checking if *keys (nested) exists in `element` (dict).
def keys_exists(element, *keys):
    if not isinstance(element, dict):
        raise AttributeError('keys_exists() expects dict as first argument.')
    if len(keys) == 0:
        raise AttributeError('keys_exists() expects at least two arguments, one given.')

    _element = element
    for key in keys:
        try:
            _element = _element[key]
        except KeyError:
            return False
    return True
node_id_count=1
with open(in_fp, "r") as f_in:
    for line in f_in:
        line=line.strip()
        line=line.split("\t")
        
        ID_A=line[0] if not line[0].startswith("WUHANREF") else line[0]
        ID_B=line[1] if not line[1].startswith("WUHANREF") else line[1]
        node_id_short_A="_".join(ID_A.split("_")[1:])
        node_id_short_B="_".join(ID_B.split("_")[1:])

        for current_ID in [ID_A, ID_B]:
            if current_ID.startswith("WUHANREF"):
                current_ID_short="_".join(current_ID.split("_")[1:])
                seq_details[current_ID_short]={}
                seq_details[current_ID_short]['gene'] = "_".join(current_ID.split("_")[1:])
                seq_details[current_ID_short]['sample'] = "WUHANREF"
                seq_details[current_ID_short]['species'] = "WUHANREF"
                seq_details[current_ID_short]['geneLength'] = "NA"
            elif current_ID.startswith("ORF10_") or current_ID.startswith("ORFE_")  or current_ID.startswith("ORF8_"):
                current_ID_short="_".join(current_ID.split("_")[1:])
                seq_details[current_ID_short]={}
                seq_details[current_ID_short]['gene'] = current_ID.split("_")[0]
                seq_details[current_ID_short]['sample'] = "_".join(current_ID.split("_")[1:])
                seq_details[current_ID_short]['species'] = current_ID.split("_")[1]
                seq_details[current_ID_short]['geneLength'] = "NA"

        sample_A=seq_details[node_id_short_A]['sample']
        sample_B=seq_details[node_id_short_B]['sample']
        
        if not sample_A == sample_B:
            query_node = line[0]
            subject_node = line[1]
            # Only record the lines where sample A is compared to sample B and the match is on the plus strand (to remove duplicate and reverse mapping)            gene_A= line[0].split("_")
            filtered_blast_results.update({line[0]:{}}) if line[0] not in filtered_blast_results else 0
            filtered_blast_results[line[0]][line[1]] = {
                                                         'pident':line[2],
                                                         'lenA':line[5], 
                                                         'lenB':line[6], 
                                                         'lenMatch':line[7], 
                                                         'lenMismatch':line[9], 
                                                         'qcovs': line[8], 
                                                         'gapOpen':line[10],
                                                         'blastScore':line[11],                                                             
                                                         'evalue':line[4]                                                                 
                                                         }

            nodes= [ID_A, ID_B]
            for node_idx, node in enumerate(nodes):
                node_name = nodes[node_idx]
                if node_name not in node_annotation:
                    node_short = "_".join(node.split("_")[1:])
                    current_gene=seq_details[node_short]['gene']
                    current_gene =  gene_name_to_wuhanref_gene_name[current_gene] if current_gene in gene_name_to_wuhanref_gene_name else current_gene
                    node_annotation[node] = {}
                    node_annotation[node]['label'] = node_name
                    node_annotation[node]['sample'] = seq_details[node_short]['sample']
                    node_annotation[node]['species'] = seq_details[node_short]['species']
                    node_annotation[node]['gene'] = current_gene
                    node_annotation[node]['geneLength'] = seq_details[node_short]['geneLength']
                    if node_name not in node_name2id:
                        node_name2id[node]=node_id_count
                        node_id2name[node_id_count]=node
                        node_id_count= node_id_count+1


###########
# Generate a network graph with all the blast score - this is to group clusters of similar sequences, thus enabling comparison across species
# Define "orphan nodes, where it is not the best matched node for another sample/species
f_out_flat = open(out_flattened_df_nodes, "w+", newline="")
f_out_flat_header=False
f_out=open(out_gml, "w+", newline="")
f_out.write("graph[\n")
# Print out node info
for current_node_id in range(1, node_id_count):
    f_out.write("  node [\n")
    f_out.write("     id %d\n"%(current_node_id))
    node_name=node_id2name[current_node_id]
    annotation_names = list(node_annotation[node_name].keys())
    annotation_names.sort()
    if not f_out_flat_header:
        f_out_flat.write("%s\n"%("\t".join(annotation_names)))
    f_out_flat.write(node_name)
    for current_annotation in annotation_names:
        f_out.write('     %s "%s"\n'%(current_annotation, node_annotation[node_name][current_annotation]))
        f_out_flat.write("\t" +  node_annotation[node_name][current_annotation])
    f_out.write("   ]\n")
f_out_flat.close()
# Print out edge info
for node_A in filtered_blast_results:
    for node_B in filtered_blast_results[node_A]:
        current_node_A_id=node_name2id[node_A]
        current_node_B_id=node_name2id[node_B]
        annotation_names = list(filtered_blast_results[node_A][node_B].keys())
        f_out.write("  edge [\n")
        f_out.write("     source %d\n"%(current_node_A_id))
        f_out.write("     target %d\n"%(current_node_B_id))
        for current_annotation in annotation_names:
            f_out.write('     %s "%s"\n'%(current_annotation, filtered_blast_results[node_A][node_B][current_annotation]))
        f_out.write("  ]\n")
f_out.write("]")
f_out.close()
###########
f_out = open(out_flattened_df_edges, "w+", newline="")

f_out.write((",".join(["node_A", 
                       "node_B", 
                       "pident", 
                       "node_A_len", 
                       "node_B_len", 
                       "alignment_len", 
                       "mismatches",
                       "qcovs",
                       "gap_openings",
                       "blast_score"])) + "\n")
# Print out results
for node_A in filtered_blast_results:
    for node_B in filtered_blast_results[node_A]:
        current_node_A_name=node_A
        current_node_B_name=node_B
        out_values=[current_node_A_name,
                    current_node_B_name,
                    filtered_blast_results[node_A][node_B]['pident'],
                    filtered_blast_results[node_A][node_B]['lenA'],
                    filtered_blast_results[node_A][node_B]['lenB'],
                    filtered_blast_results[node_A][node_B]['lenMatch'],
                    filtered_blast_results[node_A][node_B]['lenMismatch'],
                    filtered_blast_results[node_A][node_B]['qcovs'],
                    filtered_blast_results[node_A][node_B]['gapOpen'],
                    filtered_blast_results[node_A][node_B]['blastScore'],
                    filtered_blast_results[node_A][node_B]['evalue']]
        out_line = ",".join(out_values) + "\n"
        f_out.write(out_line)
f_out.close()
    

