#working_dir= "D:/git/coronahack/final_paper/tree"
#setwd(working_dir)

genome_info_dir <- "../host-data"

# install packages
if(FALSE){
	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")
	BiocManager::install("treeio")
	BiocManager::install("ggtree")
	install.packages("ape",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("caper",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("diversitree",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("geiger",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("nlme",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("OUwie",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("phangorn",repos="https://cloud.r-project.org",quiet=TRUE)
	install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
}
# library
library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(tidyr)
library(phytools)

codon_usage_clusters_fp <- "../codon_usage/input/Codon_usage_cluster.csv"
codon_usage_df <- read.csv(codon_usage_clusters_fp, stringsAsFactors=FALSE)
plot_lane_names = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O")
plot_lane_names<- tolower(plot_lane_names)
subspecies_fp <- "Genome_Host_Species.csv"
# filter out duplicated record based on length of hostspecies label
subspecies_df  <- read.csv(subspecies_fp , stringsAsFactors=FALSE)
subspecies_df$HostSpecies_len <- sapply(subspecies_df$HostSpecies, nchar)
subspecies_df <- subspecies_df[order(subspecies_df$HostSpecies_len, decreasing=TRUE),]
subspecies_df <- subspecies_df[!duplicated(subspecies_df$Label),]
subspecies_df$host_genus <- sapply(strsplit(subspecies_df$HostSpecies, " "), function(x)x[1])
# find host species label with > 10 samples
hostspecies_high_count <- aggregate(data=subspecies_df, Label~HostSpecies, length)
hostspecies_high_count <- hostspecies_high_count$HostSpecies[hostspecies_high_count$Label>10]
subspecies_df$hostspecies_highspecies <- ifelse(subspecies_df$HostSpecies %in% hostspecies_high_count, subspecies_df$HostSpecies, "")
subspecies_df$host_genus <- ifelse(subspecies_df$HostSpecies %in% hostspecies_high_count, "", subspecies_df$host_genus)
# label genera with > 10 samples
hostgenus_high_count <- aggregate(data=subspecies_df, Label~host_genus, length)
hostgenus_high_count <- hostgenus_high_count$host_genus[hostgenus_high_count$Label>10]
subspecies_df$bathost_label <- ifelse(subspecies_df$host_genus %in% hostgenus_high_count, subspecies_df$host_genus, "z Other")
subspecies_df$bathost_label <- apply(subspecies_df[,c("hostspecies_highspecies","bathost_label")], 1, function(x)ifelse(nchar(x[1])>0,x[1],x[2]))
subspecies_df$bathost_label <- factor(subspecies_df$bathost_label, levels=c(unique(subspecies_df$bathost_label)[!unique(subspecies_df$bathost_label) %in% c("Other", "")], "Other", ""))

# keep the species label for > 10 samples
subspecies_high_count <- aggregate(data=subspecies_df, Label~SubSpecies, length)
subspecies_high_count <- subspecies_high_count$SubSpecies[subspecies_high_count$Label > 10]
subspecies_df$subspecies_label <- ifelse(subspecies_df$SubSpecies %in% subspecies_high_count, subspecies_df$SubSpecies, "z Other")
subspecies_df$subspecies_label <- factor(subspecies_df$subspecies_label, levels=c(unique(subspecies_df$subspecies_label)[!unique(subspecies_df$subspecies_label) %in% c("Other", "")], "Other", ""))

# read in tree
nex <- "RAxML_bestTre"
nex_small <- "RAxML_bestTree.small_bat"
tree <- read.tree(nex)
tree <- unroot(tree)
tree <- midpoint.root(tree)
tree <- ladderize(tree)

tree_small <- read.tree(nex_small)
tree_small <- unroot(tree_small)
tree_small <- midpoint.root(tree_small)
tree_small <- ladderize(tree_small)
#tree <- root(tree, outgroup = "KJ473808")
#tree <- root(tree, outgroup = "GU190215")
# the sample names for the trees
all_tree_samples <- tree$tip.label


# Organise background information
genome_metrics_fp <- paste(genome_info_dir, list.files(genome_info_dir,pattern = "\\genome_metrics.csv$", recursive = T), sep="/")
genome_metrics_names <- sapply(genome_metrics_fp , FUN=function(x){out <- basename(x); return(strsplit(out, "_")[[1]][1])})
genome_metrics <- lapply(genome_metrics_fp, read.csv, stringsAsFactors=FALSE, header=FALSE)
genome_metrics <- lapply(genome_metrics, FUN=function(x){out <- x; out$V18=NULL; colnames(out) <- out[1,]; out <- out[2:nrow(out),]; return(out)})
names(genome_metrics) <- genome_metrics_names
genome_metrics <- lapply(genome_metrics_names, FUN=function(x){out <- genome_metrics[[x]];out$dataset_name <- x; return(out)})
genome_metrics_df <- do.call(rbind, genome_metrics)
write.csv(genome_metrics_df, "../gene_summary/all_genome_metrics.csv", row.names=FALSE)
row.names(genome_metrics_df) <- genome_metrics_df$Genome
# There are some mislabelling. Tidy those
genome_metrics_df$Family <- ifelse(genome_metrics_df$Family %in% c("N/A", "-N/A-"), "", genome_metrics_df$Family)
genome_metrics_df<- genome_metrics_df[!(row.names(genome_metrics_df) %in% "    `"),] 
all_tree_samples[!(all_tree_samples %in% genome_metrics_df$Genome)] #find the name of missing annotation

# genome/genes summary
genome_metrics_list <- split(genome_metrics_df, genome_metrics_df$dataset_name)



# Organise annotation_information
heatmapData <- genome_metrics_df[,c('dataset_name','Family')]
heatmapData <- heatmapData[heatmapData$dataset_name %in% c("bat", "pangolin", "wuhan"),]
heatmapData <- rbind(heatmapData, c("SARS-CoV2-reference","Betacoronavirus"))
row.names(heatmapData)[nrow(heatmapData)] <- "MN908947.3"
heatmapData$dataset_name <- gsub("wuhan", "SARS-CoV2", heatmapData$dataset_name)


# Read in VEP results
vep_bat <- read.csv("combined_variants.vep.bcsq.annot_bat.txt",stringsAsFactors=FALSE)
vep_pangolin <- read.csv("combined_variants.vep.bcsq.annot_pangolin.txt", stringsAsFactors=FALSE)
variants_interest <- c(
						"240PE>240P",
						"68S>68SD",
						"68S>68SE",
						"68S>68SQ",
						"3DS>3D",
						"30DY>30D",
						"93V>93VY",
						"93V>93VH",
						"93V>93VQ",
						"7Q>7QP",
						"7Q>7QS",
						"26Y>26*"						
						)# the variants found in the table 
# these are for grouping variants by gene instead of position
variants_interest_2 <- c("927PD>927P" , 
						"1227Q>1227A" , 
						"1227Q>1200QA" , 
						"3164R>3164RR",
						"3573KR>3574K",
						"3576I>3575IV",
						"6P>6PQ",
						"238GQ>239G",
						"385RQ>385R"
							)

vep <- rbind(vep_bat, vep_pangolin)
vep <- unique(vep[vep$amino_acid_change %in% variants_interest ,])
vep <- split(vep, vep$POS)

vep_unique <- rbind(vep_bat, vep_pangolin)
vep_unique <- unique(vep_unique[,c("gene", "CHROM","POS", "REF", "ALT", "consequence" )])
consequence_unique <- unique(vep_unique$consequence)
consequence_unique <- consequence_unique[order(consequence_unique)]
consequence_kept <- c("inframe_deletion", "inframe_insertion", "missense", "stop_gained", "stop_lost")
vep_unique <- unique(vep_unique[vep_unique$consequence %in% consequence_kept, ])



extract_vep_samples <- function(f_df){
	f_df$gene_consequence_aa <- paste(f_df$gene,f_df$consequence, f_df$amino_acid_change, sep="-")
	f_df_sample_list <- split(f_df$VariantSamples, f_df$gene_consequence_aa)
	f_out <- lapply(f_df_sample_list, FUN=function(x)strsplit(x, "\\|")[[1]])
	f_out <- lapply(names(f_out), FUN=function(x){out <- data.frame(Genome=f_out[[x]], change=x, stringsAsFactors=FALSE)})
	if(is.list(f_out)){
		f_out <- do.call(rbind, f_out)
	}
	return(f_out)
}
vep_selected <- lapply(vep, extract_vep_samples)
vep_selected <- do.call(rbind, vep_selected)
vep_selected$change_position <- sapply(strsplit(vep_selected$change, ">"), FUN=function(x)x[1]) 
#vep_selected$gene <- sapply(strsplit(vep_selected$change, "-"), FUN=function(x)x[1]) 
vep_selected <- spread(vep_selected, change_position, change)

vep_2 <- rbind(vep_bat, vep_pangolin)
vep_2 <- unique(vep_2[vep_2$amino_acid_change %in% variants_interest_2 ,])
vep_2 <- split(vep_2, vep_2$POS)
vep_selected_2 <- lapply(vep_2, extract_vep_samples)
vep_selected_2 <- do.call(rbind, vep_selected_2)
vep_selected_2$gene <- sapply(strsplit(vep_selected_2$change, "-"), FUN=function(x)x[1]) 
vep_selected_2$change_position <- sapply(strsplit(vep_selected_2$change, ">"), FUN=function(x)x[1]) 
## two variants on the same gene for the same 2 samples, manually merge them into one entry
#merging_entries <- c("ORF1ab-*inframe_insertion-1227Q>1200QA", "ORF1ab-inframe_deletion-927PD>927P")
#output_entery <- "ORF1ab-inframe_insertion-1227Q>1200QA-inframe_deletion-927PD>927P"
#vep_selected_2_sub <- vep_selected_2[vep_selected_2$change %in% merging_entries,]
#vep_selected_2_sub$change <- output_entery

# remove this specific variant ORF1ab-*inframe_insertion-1227Q>1200QA
vep_selected_2 <- vep_selected_2[!(vep_selected_2$change %in% "ORF1ab-*inframe_insertion-1227Q>1200QA"), ]
vep_selected_2$change_position = NULL
vep_selected_2 <- spread(vep_selected_2, gene, change)
vep_selected <- merge(vep_selected, vep_selected_2, by="Genome", all.x=TRUE,all.y=TRUE)



# combine all ring annotation 
heatmapData_2 <- merge(heatmapData, subspecies_df[,c('Label','subspecies_label','bathost_label')], by.x=0, by.y="Label", all.x=TRUE)
heatmapData_2 <- heatmapData_2[,c("Row.names", "dataset_name", "bathost_label", "Family", "subspecies_label")]
heatmapData_2 <- merge(heatmapData_2, codon_usage_df[,c("Genome","cluster")], by.x="Row.names", by.y="Genome", all.x=TRUE)
vep_multi_var_same_pos <- c("Genome", "E-inframe_insertion-68S", "N-inframe_insertion-7Q", "ORF7a-inframe_insertion-93V")
vep_selected <- vep_selected[,c(vep_multi_var_same_pos, colnames(vep_selected)[!(colnames(vep_selected) %in% vep_multi_var_same_pos)])]
heatmapData_2 <- merge(heatmapData_2, vep_selected, by.x="Row.names", by.y="Genome", all.x=TRUE)
row.names(heatmapData_2) <- heatmapData_2$Row.names
heatmapData_2$Row.names=NULL
all_unique <- do.call(c, sapply(heatmapData_2, function(x) (unique(x[!is.na(x)]))))
all_unique <- all_unique[order(all_unique)]
heatmapData_2_original <- heatmapData_2
heatmapData_2 <- heatmapData_2[,!(colnames(heatmapData_2) %in% c("Family", "subspecies_label"))]

for(i in 1:ncol(heatmapData_2)){
	current_col_name=colnames(heatmapData_2)[i]
	current_lane_name=plot_lane_names[i]
	NA_val <- is.na(heatmapData_2[[current_col_name]])
	space_val <- heatmapData_2[[current_col_name]] == ""
	empty_val <- NA_val | space_val
	heatmapData_2[[current_col_name]] <- paste(current_lane_name, heatmapData_2[[current_col_name]], sep='  ')
	heatmapData_2[[current_col_name]][empty_val] = ""
	heatmapData_2[[current_col_name]] <- gsub("_"," ", heatmapData_2[[current_col_name]])

}

genome_metrics_df_plot <- merge(genome_metrics_df, heatmapData_2_original[,c("bathost_label","subspecies_label")], by.x="Genome", by.y=0, all.x=TRUE)
genome_metrics_df_plot$subspecies_label <- as.character(genome_metrics_df_plot$subspecies_label)
genome_metrics_df_plot$subspecies_label <- ifelse(genome_metrics_df_plot$subspecies_label %in% c("Sarbecovirus", "Nobecovirus", "Merbecovirus"), genome_metrics_df_plot$subspecies_label,NA )
genome_metrics_df_plot$subspecies_label <- ifelse(genome_metrics_df_plot$dataset_name %in% c("wuhan", "charite"), "Sarbecovirus", genome_metrics_df_plot$subspecies_label)
genome_metrics_df_plot$Genome2 <- genome_metrics_df_plot$Genome
# Draw tree (all)

ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
ggsave("tree_node_number.png", height=1200,width=1200,units='mm') # use this to find nodes where swapping the orientation might look better aethetically 
p <- ggtree(tree,  layout="circular") %<+% genome_metrics_df_plot 
p <- ggtree(tree) %<+% genome_metrics_df_plot  
# swap these two nodes
p <- rotate(p, 483)


#p <- ggtree(tree) %<+% genome_metrics_df_plot  +
#    aes(color=subspecies_label) 
#ggsave("test_big.png", height=1200,width=400,units='mm')
#p <- ggtree(tree,  layout="circular") %<+% genome_metrics_df_plot +
#	geom_tippoint(aes(color=dataset_name)) 
#	geom_tippoint(aes(color=dataset_name)) +
#	geom_tiplab2() 
# Add annotation
heatmap.colours <- c("#33a02c", #"bat" 1_1
					"#ff7f00",	#"pangolin" 1_2
					"#6a3d9a",	# "SARS-CoV2" 1_3
					"#a6cee3",	#"SARS-CoV2-reference" 1_4

					"#33a02c", # host species -- Manis
					"#B0E2FF", #Rhinolophus
					"#9370DB", #Rhinolophus ferr
					"#104E8B", #Rhinolophus sins
					"#fb9a99", #Scotophilus
					"#AAAAAA", #Other

#					"#33a02c",  #"Alphacoronavirus" 2_1
#					"#ff7f00", 	#"Betacoronavirus" 2_2

#					"#33a02c", #sub speceis -- Decacovirus
#					"#ff7f00", #Merbecovirus
#					"#6a3d9a", #Minunacovirus
#					"#a6cee3", #Nobecovirus 
#					"#1f78b4", #Rhinacovirus
#					"#e31a1c", #Sarbecovirus
#					"#AAAAAA", #Other
					
                    "#1f78b4",	#"codon_usage_cluster1
                    "#ff0000",	#"codon_usage_cluster2
                    "#b2df8a",	#"codon_usage_cluster3

					"#843b62", 	#"E-inframe_insertion-68S>68SD"  3_1
					"#fb9a99", 	#"E-inframe_insertion-68S>68SE" 3_2
					"#e31a1c", 	#"E-inframe_insertion-68S>68SQ"  4_3
					"#33a02c", 	#"N-inframe_insertion-7Q>7QP" 6_1
					"#b2df8a",	#"N-inframe_insertion-7Q>7QS" 6_2
                    "#ffc300",	#"ORF7a-inframe_insertion-93V>93VY" 9_1
                    "#900c3f",	#"ORF7a-inframe_insertion-93V>93VH" 9_2
                    "#ff5733",	#"ORF7a-inframe_insertion-93V>93VQ" 9_3
                    "#b15928", 	#"M-inframe_deletion-3DS>3D" 5_1
					"#48C9B0",	#"ORF10-stop_gained-26Y>26*" 7_1
					"#ffff99",	#"ORF3a-inframe_deletion-240PE>240P" 8_1
					"#b2df8a", 	#"ORF6-inframe_deletion-30DY>30D" 6_1
					"#843b62", 			#"N-inframe_deletion-385RQ>385R"
					"#fb9a99", 			#"N-inframe_deletion-238GQ>239G"
					"#e31a1c", 			#"N-inframe_insertion-6P>6PQ"
					"#fee08b", 			#"ORF1ab-inframe_deletion-3573KR>3574K"
					"#d53e4f", 			#"ORF1ab-inframe_deletion-927PD>927P" 
					"#99d594", 			#"ORF1ab-inframe_insertion-3164R>3164RR"
					"#3288bd", 			#"ORF1ab-inframe_insertion-3576I>3575IV"				
					"#ffffff"
					)

p2 <- gheatmap(p, heatmapData_2, offset = 0, color=NULL, 
         colnames_position="top",  width=0.6,
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=0.5)  +
		 scale_fill_manual(values=heatmap.colours) 

ggsave("all_tree.png", width=450,height=300, unit='mm')
ggsave("all_tree.pdf", width=11, height=8.5)


heatmap.colours2 <- c("#33a02c", #"bat" 1_1
					"#ff7f00",	#"pangolin" 1_2
					"#6a3d9a",	# "SARS-CoV2" 1_3
					"#a6cee3",	#"SARS-CoV2-reference" 1_4

#					"#33a02c", # host species -- Manis
#					"#B0E2FF", #Rhinolophus
#					"#9370DB", #Rhinolophus ferr
#					"#104E8B", #Rhinolophus sins
#					"#fb9a99", #Scotophilus
#					"#AAAAAA", #Other

#					"#33a02c",  #"Alphacoronavirus" 2_1
#					"#ff7f00", 	#"Betacoronavirus" 2_2

#					"#33a02c", #sub speceis -- Decacovirus
#					"#ff7f00", #Merbecovirus
#					"#6a3d9a", #Minunacovirus
#					"#a6cee3", #Nobecovirus 
#					"#1f78b4", #Rhinacovirus
#					"#e31a1c", #Sarbecovirus
#					"#AAAAAA", #Other
					
                    "#1f78b4",	#"codon_usage_cluster1
                    "#ff0000",	#"codon_usage_cluster2
                    "#b2df8a",	#"codon_usage_cluster3

					"#843b62", 	#"E-inframe_insertion-68S>68SD"  3_1
					"#fb9a99", 	#"E-inframe_insertion-68S>68SE" 3_2
					"#e31a1c", 	#"E-inframe_insertion-68S>68SQ"  4_3
					"#33a02c", 	#"N-inframe_insertion-7Q>7QP" 6_1
					"#b2df8a",	#"N-inframe_insertion-7Q>7QS" 6_2
                    "#ffc300",	#"ORF7a-inframe_insertion-93V>93VY" 9_1
                    "#900c3f",	#"ORF7a-inframe_insertion-93V>93VH" 9_2
                    "#ff5733",	#"ORF7a-inframe_insertion-93V>93VQ" 9_3
                    "#b15928", 	#"M-inframe_deletion-3DS>3D" 5_1
					"#48C9B0",	#"ORF10-stop_gained-26Y>26*" 7_1
					"#ffff99",	#"ORF3a-inframe_deletion-240PE>240P" 8_1
					"#b2df8a", 	#"ORF6-inframe_deletion-30DY>30D" 6_1
					"#843b62", 			#"N-inframe_deletion-385RQ>385R"
					"#fb9a99", 			#"N-inframe_deletion-238GQ>239G"
					"#e31a1c", 			#"N-inframe_insertion-6P>6PQ"
					"#fee08b", 			#"ORF1ab-inframe_deletion-3573KR>3574K"
					"#d53e4f", 			#"ORF1ab-inframe_deletion-927PD>927P" 
					"#99d594", 			#"ORF1ab-inframe_insertion-3164R>3164RR"
					"#3288bd", 			#"ORF1ab-inframe_insertion-3576I>3575IV"				
					"#ffffff"
					)

heatmapData_3 <- heatmapData_2[,colnames(heatmapData_2)[!(colnames(heatmapData_2) %in% c("bathost_label"))]]

p3 <- gheatmap(p, heatmapData_3, offset = 0.6, color=NULL, 
         colnames_position="top",  width=0.2,
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=0.5)  +
		 scale_fill_manual(values=heatmap.colours2) 
p3 <- viewClade(p3+geom_tiplab(), node=381)


ggsave("all_tree_clade.png", width=650,height=400, unit='mm')
ggsave("all_tree_clade.pdf", width=380,height=400, unit='mm')





