library(tidyr)
library(ca)
library(ggplot2)
library(gridExtra)
library(MCL)

#working_dir="D:/git/coronahack/final_paper/codon_usage"
#setwd(working_dir)

#### Settings
out_rscu_fp <- "input/rscu_per_gene.csv"
in_fp <- "input/codon_usage_best_gene.csv" # this is produced by script/find_best_blast_hit.py
input_dir="input/codon_count"
inupt_gene_count <- "output/summary/sample_gene_count.csv"
output_dir="output/pca"

dir.create(output_dir, showWarnings=FALSE)
genome_info_dir <- "../host-data"
dir.create(dirname(input_dir), showWarnings=FALSE)
dir.create(output_dir, showWarnings=FALSE)


#### Organise background information
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



##### generate aa tables
acids<-c("Isoleucine","Leucine","Valine","Phenylalanine","Methionine","Cysteine","Alanine","Glycine","Proline","Threonine","Serine",
	 "Tyrosine","Tryptophan","Glutamine","Asparagine","Histidine","Glutamic acid","Aspartic acid","Lysine","Arginine","Stop codons")
slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","Stop")
codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
	 "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
	 "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")
ignore_codons <- c("TAA", "TGA", "TAG", "AUG", "TGG") # ignore these start, stop and aa that code from one codon
aa_table <- data.frame(aa=acids,slc=slc,codon=codon,stringsAsFactors=FALSE)
aa_table <- split(aa_table, aa_table$aa)
aa_table_per_codon <-  do.call(rbind,lapply(aa_table, FUN=function(x)data.frame(codon=strsplit(x$codon, ", ")[[1]],aa=x$aa[1], stringsAsFactors=FALSE)))
aa_table_per_codon <- aa_table_per_codon[!(aa_table_per_codon$codon %in% ignore_codons),]
aa_table_per_codon_max_count <- aggregate(data=aa_table_per_codon, codon~aa, length)
colnames(aa_table_per_codon_max_count) <- c("aa","max_num_codon")
fill_col=c("#33a02c","#ff7f00", "#6a3d9a","#1f78b4", "#BDBDBD")
fill_col2 <- rep(c("#ffffff", "#D3D3D3", "#7D7D7D", "#000000"), 64/4)
fill_order = c("bat", "pangolin", "wuhan", "charite", "codon") 
fill_order2 = c("BAT", "PANGOLIN", "WUHAN", "CHAIRITE", "Codon") 

#### functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
plot_codon_bias <- function(f_plot_df, f_name){
	f_plot_df$dataset_name <- factor(full_df$dataset_name, levels=c("Bat", "Pangolin", "Charite", "Wuhan"))
	p <- ggplot(f_plot_df, aes(x=aa, y=codon_count, fill=codon)) + 
		geom_bar(position="fill", stat="identity") + 
		scale_fill_manual(values=fill_col2) + 
		theme(legend.position = "none") +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		facet_grid(dataset_name ~ .) 
	ggsave(paste(output_dir, "/", f_name, ".png", sep=""))
	return(p)
}

calculate_codon_ratio_rscu <- function(y){
	f_df_split <- split(y, y$aa)
	f_rscu <- lapply(f_df_split, function(z){z_out <- z; z_out$rscu <- z_out$codon_count/sum(z_out$codon_count)*length(aa_table_per_codon$aa[aa_table_per_codon$aa %in% z$aa[1]])  ; return(z_out)})
	f_rscu <- lapply(f_rscu, function(z){z_out <- z; z_out$codon_ratio <- z_out$codon_count/sum(z_out$codon_count)  ; return(z_out)})
	f_rscu <- do.call(rbind, f_rscu)
	# make sure that there is a reading for each codon (i.e. fill the missing amino acids with zero)
	missing_codon <- aa_table_per_codon[!(aa_table_per_codon$codon %in% f_rscu$codon),]
	if(nrow(missing_codon)>0){
		missing_codon$Genome=f_rscu$Genome[1]
		missing_codon$gene_sample_name=f_rscu$gene_sample_name[1]
		missing_codon$gene=f_rscu$gene[1]
		missing_codon$dataset_name =f_rscu$dataset_name[1]
		missing_codon$codon_count =0
		missing_codon$codon_ratio =0
		missing_codon$codon_rscu =0
		missing_codon <- missing_codon[,colnames(f_rscu)]
		f_rscu <- rbind(f_rscu, missing_codon)
	}
	f_rscu <- f_rscu[f_rscu$codon %in% aa_table_per_codon$codon,]
	return(f_rscu)
}
bad_samples = "" #c('KC881005','KC881005') # need to check
plot_genes = c( 'more_than_half_all_species',  'ORF1ab','S','ORF3a','E','ORF6','ORF7a','ORF7b','N', 'ORF10', 'ORF8','M', 'ORF1a')
all_genes = c('ORF1ab','S','ORF3a','E','ORF6','ORF7a','ORF7b','N', 'ORF10', 'ORF8','M', 'ORF1a')


### loop through to calcuate and plot RSCU for each gene/gene group
for (current_gene in plot_genes){
	print(current_gene)
	genome_collpase = current_gene # change this to no_orf8/all/core
	##### gene selection
	if(genome_collpase %in% "no_orf8"){
		keep_genes = c('ORF1ab','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','N', 'ORF10')
	} else if (genome_collpase %in% "no_orf86"){
		keep_genes = c('ORF1ab','S','ORF3a','E','M','ORF7a','ORF7b','N', 'ORF10')
	} else if (genome_collpase %in% "all"){
		keep_genes = c('ORF1ab','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','N', 'ORF10', 'ORF8')
	} else if (genome_collpase %in% "core"){
		keep_genes = c('E','M','N', 'S')
	} else if (genome_collpase %in% "more_than_half_all_species"){
		keep_genes = c('E','N', 'S', 'ORF1a', 'ORF3a', 'ORF10') # removed M
		#keep_genes = c('E','M','N', 'S', 'ORF1ab', 'ORF3a', 'ORF10')
	} else if (genome_collpase %in% "core_1ab"){
		keep_genes = c('E','M','N', 'S', 'ORF1ab')
	} else {
		keep_genes = genome_collpase
	}



	#### These are done by aggregating "input/codon_usage_best_gene.csv" produced by find_best_blast_hit.py
	raw_data <- read.csv(in_fp, stringsAsFactors=FALSE)
	raw_data$gene<- sapply(strsplit(as.character(raw_data$gene_sample_name), "_"), FUN=function(x)x[1])
	raw_data$Genome <- sapply(strsplit(as.character(raw_data$gene_sample_name), "_"), FUN=function(x)paste(x[2:length(x)],collapse="_"))
	raw_data <- raw_data[raw_data$gene %in% keep_genes,]

	# count number of genes per genome, keep only genomes with max count if it is not the "all" selection
	raw_data_Genome_count <- aggregate(data=raw_data, gene~Genome, length)
	colnames(raw_data_Genome_count) <- c("Genome", "gene_count")
	raw_data <- merge(raw_data, raw_data_Genome_count, by="Genome")
	raw_data <- raw_data[raw_data$gene_count == max(raw_data$gene_count),]
	raw_data <- gather(raw_data, codon, codon_count, AAA:YTT)
	full_df <- merge(raw_data, aa_table_per_codon, by="codon")
	full_df <- merge(full_df, genome_metrics_df[,c("Genome","dataset_name")], by="Genome")
	full_df <- full_df[!(full_df$Genome %in% bad_samples), ]
	# merge all the data together so is for the full genome (rather than for each gene)
	full_df_agg <- aggregate(data=full_df, codon_count~codon+dataset_name + aa + Genome, sum)

	# calculate rscu
	full_df_agg_list <- split(full_df_agg, full_df_agg$Genome)
	full_df_agg_list <- lapply(full_df_agg_list, calculate_codon_ratio_rscu)
	full_df_agg <- do.call(rbind, full_df_agg_list)
	full_df_agg_dataset <- aggregate(data= full_df_agg, rscu~codon + dataset_name + aa, mean)
	full_df_agg_view <- spread(full_df_agg_dataset, dataset_name, rscu)
	full_df_agg_view <- full_df_agg_view[order(full_df_agg_view$aa, full_df_agg_view$codon),]
	

	
	full_df_ca <- spread(full_df_agg[,c("Genome", "rscu", "codon")], Genome, rscu)
	sample_info <- unique(full_df_agg[,c("Genome", "dataset_name")])
	row.names(sample_info) <- sample_info$Genome
	# plot Correspondence Analysis for rscu across all core genes
	ca_input <- (full_df_ca[,2:ncol(full_df_ca)])
	row.names(ca_input) <- full_df_ca$codon
	
	
	# remove columns with no standard deviation
	check_zero_sd  <- apply(ca_input,1,sd)
	pca_mx <- t(ca_input[check_zero_sd>0,])
	pca_mx[is.na(pca_mx)]<-0
	colnames(pca_mx) <- rownames(ca_input[check_zero_sd>0,])
	rownames(pca_mx) <- colnames(ca_input)
	
	pca <- prcomp(pca_mx)
 	eigs <- pca$sdev^2
	PC1_cont <-  format((eigs[1] / sum(eigs) *100), digits=2, nsmall=2)
	PC2_cont <-  format((eigs[2] / sum(eigs) *100), digits=2, nsmall=2)
	
	pca_plot_df <- as.data.frame(pca$x[,1:2])
	pca_plot_df <- merge(pca_plot_df, sample_info, by=0)
	pca_plot_df$Row.names=NULL
	pca_plot_df$dataset_name <- factor(pca_plot_df$dataset_name, levels=fill_order)
	# label the three odd bat genomes
	pca_plot_df$label <- ifelse(pca_plot_df$Genome %in%  c("MN996532","MG772934","MG772933"), pca_plot_df$Genome , "")
	if(genome_collpase %in% all_genes){
		plot_title_label=genome_collpase
	}else{
		plot_title_label=""
	}
#	if(paste0(genome_collpase, collapse="") == "ORF1a"){
#		ORF1a_outlier = "MG762674"
#		pca_plot_df<- pca_plot_df[!(pca_plot_df$Genome %in% ORF1a_outlier), ]
#	}
	# label the total number of genes for each dataset at the top left corner
	# count the total number of genomes per dataset
	genome_count = aggregate(data=unique(pca_plot_df[,c("Genome", "dataset_name")]), Genome~dataset_name, length)
	count_bat=genome_count$Genome[genome_count$dataset_name %in% "bat"]
	count_charite=genome_count$Genome[genome_count$dataset_name %in% "charite"]
	count_wuhan=genome_count$Genome[genome_count$dataset_name %in% "wuhan"]
	count_pangolin=genome_count$Genome[genome_count$dataset_name %in% "pangolin"]
	plot_count_label = paste("\n\n\n\n\n b:", count_bat, "\n", 
							" p:", count_pangolin, "\n", 
							" c:", count_charite, "\n", 
							" w:", count_wuhan, "\n", sep="")
	ggplot(pca_plot_df, aes(x=PC1, y=PC2, col = dataset_name, label=label)) + 
			geom_point(alpha=0.7, size = 2) +
			theme_bw() +
			geom_text(size=3,hjust = 0, nudge_x = 0.1)+
			scale_colour_manual(values=fill_col) + 
			xlab(paste0("PC1 (",PC1_cont, "%)")) +
			ylab(paste0("PC2 (",PC2_cont, "%)")) +
			annotate("text",  x=-Inf, y = Inf, label = plot_count_label, size=3.5, hjust = 0)	+
			theme(legend.position = "none")
	
	ggsave(paste0(output_dir, "/", genome_collpase, ".png"),width=80,height=80, unit='mm')

	# get kmeans clusters if the plot genes are more_than_half_all_species
	if( "more_than_half_all_species" %in% genome_collpase){
		# takes just the bat (minus RaTG13 - MN996532, MG772934 and MG772933) and perform kmeans clustering
		multiple_genes <- full_df_agg
		bat_genomes_to_kmean_cluster <- sample_info$Genome[sample_info$dataset_name %in% "bat"]
		bat_genomes_to_kmean_cluster <- bat_genomes_to_kmean_cluster[!(bat_genomes_to_kmean_cluster %in% c("MN996532","MG772934","MG772933"))]
		# k mean cluster
		set.seed(50) # kmeans clustering gives slightly different results/clustering each time, so I have used setseed to make sure that it generates the same result each time
		kmeans_cluster <- kmeans(pca$x[bat_genomes_to_kmean_cluster,], centers=3)
		kmeans_cluster <- data.frame(cluster = paste("Codon_usage_cluster", kmeans_cluster$cluster, sep="_"), Genome = names(kmeans_cluster$cluster), stringsAsFactors=FALSE)
		write.csv(kmeans_cluster, "input/codon_usage_cluster.csv")
		kmeans_cluster_pca_plot_df <- merge(pca_plot_df, kmeans_cluster, by="Genome", all.x=TRUE)
		kmeans_cluster_pca_plot_df$cluster <- ifelse(is.na(kmeans_cluster_pca_plot_df$cluster), "NA", kmeans_cluster_pca_plot_df$cluster)

		# label the number of genomes in cluster 1 and 2
		cluster_1_genome_number=length(kmeans_cluster_pca_plot_df$cluster[kmeans_cluster_pca_plot_df$cluster %in% "Codon_usage_cluster_1"])
		cluster_2_genome_number=length(kmeans_cluster_pca_plot_df$cluster[kmeans_cluster_pca_plot_df$cluster %in% "Codon_usage_cluster_2"])
		cluster_3_genome_number=length(kmeans_cluster_pca_plot_df$cluster[kmeans_cluster_pca_plot_df$cluster %in% "Codon_usage_cluster_3"])
		plot_count_label = paste0(plot_count_label, "cl1:",  cluster_1_genome_number, "\n")		
		plot_count_label = paste0(plot_count_label, "cl2:",  cluster_2_genome_number, "\n")		
		plot_count_label = paste0(plot_count_label, "cl3:",  cluster_3_genome_number)	
		kmeans_cluster_pca_plot_df$label <- ifelse(kmeans_cluster_pca_plot_df$Genome %in%  c("KY352407"), kmeans_cluster_pca_plot_df$Genome , "")
		# manual fill the colours
		cluster_fill_col <- c("#1f78b4", "#ff0000","#b2df8a", "#dddbeb")
		ggplot(kmeans_cluster_pca_plot_df, aes(x=PC1, y=PC2, col = cluster, label=label )) + 
				geom_point(alpha=1, size = 2) +
				theme_bw() +
				geom_text(size=3,hjust = 0, nudge_x = 0.1)+
				scale_colour_manual(values=cluster_fill_col) + 
				xlab(paste0("PC1 (",PC1_cont, "%)")) +
				ylab(paste0("PC2 (",PC2_cont, "%)")) +
				annotate("text",  x=-Inf, y = Inf, label = plot_count_label, size=3.5, hjust = 0)	+
				theme(legend.position = "none")
		ggsave(paste0(output_dir, "/", genome_collpase, "_cluster.png"),width=80,height=80, unit='mm')
#		ggsave(paste0(output_dir, "/", genome_collpase, "_cluster.pdf"))
		ggplot(kmeans_cluster_pca_plot_df, aes(x=PC1, y=PC2, col = cluster, label=Genome )) + 
				geom_point(alpha=1, size = 2) +
				theme_bw() +
				geom_text() +
				scale_colour_manual(values=cluster_fill_col) + 
				xlab(paste0("PC1 (",PC1_cont, "%)")) +
				ylab(paste0("PC2 (",PC2_cont, "%)")) +
				annotate("text",  x=-Inf, y = Inf, label = plot_count_label, size=3.5, hjust = 0)	+
				theme(legend.position = "none")
		ggsave(paste0(output_dir, "/", genome_collpase, "_cluster_with_label.png"),width=80,height=80, unit='mm')
	} else {
		kmeans_cluster_pca_plot_df <- merge(pca_plot_df, kmeans_cluster, by="Genome", all.x=TRUE)
		# because not all genomes were in the original kmeans clustering, those that were not are unlabelled bat cov
		kmeans_cluster_pca_plot_df$cluster <- ifelse(is.na(kmeans_cluster_pca_plot_df$cluster), "NA", kmeans_cluster_pca_plot_df$cluster)
		kmeans_cluster_pca_plot_df$cluster <- ifelse(((kmeans_cluster_pca_plot_df$cluster %in% "NA") & (kmeans_cluster_pca_plot_df$dataset_name %in% "bat") & !(kmeans_cluster_pca_plot_df$Genome %in% c("MN996532","MG772934","MG772933"))), "unclassified_bat", kmeans_cluster_pca_plot_df$cluster)
		# because not all cluster may have members present, loop through to find out what fill colour/factor levels to use
		cluster_color=c("#1f78b4" , "#ff0000",  "#b2df8a",  "#dddbeb", "#f7f4be" )
		cluster_label=c("Codon_usage_cluster_1" , "Codon_usage_cluster_2",  "Codon_usage_cluster_3",  "NA", "unclassified_bat" )
		current_cluster_levels = vector()
		current_cluster_col = vector()
		for(idx in 1:length(cluster_label)){
			if(cluster_label[idx] %in% kmeans_cluster_pca_plot_df$cluster){
				current_cluster_levels <- c(current_cluster_levels, cluster_label[idx])
				current_cluster_col <- c(current_cluster_col, cluster_color[idx])
			}
		}
		kmeans_cluster_pca_plot_df$cluster <- factor(kmeans_cluster_pca_plot_df$cluster, levels = current_cluster_levels )
		kmeans_cluster_pca_plot_df$label <- ifelse(kmeans_cluster_pca_plot_df$Genome %in%  c("KY352407"), kmeans_cluster_pca_plot_df$Genome , "")
		# label the number of genomes in cluster 1 and 2
		cluster_1_genome_number=length(kmeans_cluster_pca_plot_df$cluster[kmeans_cluster_pca_plot_df$cluster %in% "Codon_usage_cluster_1"])
		cluster_2_genome_number=length(kmeans_cluster_pca_plot_df$cluster[kmeans_cluster_pca_plot_df$cluster %in% "Codon_usage_cluster_2"])
		cluster_3_genome_number=length(kmeans_cluster_pca_plot_df$cluster[kmeans_cluster_pca_plot_df$cluster %in% "Codon_usage_cluster_3"])
		plot_count_label = paste0("\n\n", plot_count_label) 
		plot_count_label = paste0(plot_count_label, "cl1:",  cluster_1_genome_number, "\n")		
		plot_count_label = paste0(plot_count_label, "cl2:",  cluster_2_genome_number, "\n")		
		plot_count_label = paste0(plot_count_label, "cl3:",  cluster_3_genome_number)	
		kmeans_cluster_pca_plot_df <- kmeans_cluster_pca_plot_df[order(kmeans_cluster_pca_plot_df$cluster),]
		ggplot(kmeans_cluster_pca_plot_df, aes(x=PC1, y=PC2, col = cluster, label=label )) + 
				geom_point(alpha=0.5, size = 2) +
				theme_bw() +
				scale_colour_manual(values=current_cluster_col) + 
				geom_text(size=3,hjust = 0, nudge_x = 0.1)+
				xlab(paste0("PC1 (",PC1_cont, "%)")) +
				ylab(paste0("PC2 (",PC2_cont, "%)")) +
				annotate("text",  x=-Inf, y = Inf, label = plot_count_label, size=3.5, hjust = 0)	+
				theme(legend.position = "none")
		ggsave(paste0(output_dir, "/", genome_collpase, "_cluster.png"),width=80,height=80, unit='mm')
#		ggsave(paste0(output_dir, "/", genome_collpase, "_cluster.pdf"))
	}
}

# Plot a heatmap of the codon usage
#### These are done by aggregating "input/codon_usage_best_gene.csv" produced by find_best_blast_hit.py


# this is the rscu from multiple genes 
multiple_genes$gene <- "multiple genes*"
raw_data <- read.csv(in_fp, stringsAsFactors=FALSE)
raw_data$gene<- sapply(strsplit(as.character(raw_data$gene_sample_name), "_"), FUN=function(x)x[1])
raw_data$Genome <- sapply(strsplit(as.character(raw_data$gene_sample_name), "_"), FUN=function(x)paste(x[2:length(x)],collapse="_"))
raw_data <- gather(raw_data, codon, codon_count, AAA:YTT)
full_df <- merge(raw_data, aa_table_per_codon, by="codon")
full_df <- merge(full_df, genome_metrics_df[,c("Genome","dataset_name")], by="Genome")
full_df <- full_df[!(full_df$Genome %in% bad_samples), ]
by_gene_count <- split(full_df, full_df$gene)
by_gene_count <- lapply(by_gene_count, function(x)split(x, x$Genome))
by_gene_count <- lapply(by_gene_count, function(x)lapply(x, function(z)calculate_codon_ratio_rscu(z)))
by_gene_count <- do.call(rbind,lapply(by_gene_count, function(x)do.call(rbind,x)))
by_gene_count$gene_sample_name = NULL
by_gene_count <- rbind(by_gene_count, multiple_genes) # add in the values for multiple genes
by_gene_count <- merge(by_gene_count, aa_table_per_codon_max_count, by="aa")
by_gene_count <- by_gene_count[by_gene_count$max_num_codon>1,]
by_gene_plot_count <- aggregate(data=unique(by_gene_count[,c("Genome", "gene", "dataset_name")]), Genome~gene + dataset_name, length)
colnames(by_gene_plot_count) <- c("gene", "dataset_name", "genome_count")
by_gene_plot <- aggregate(data=by_gene_count, codon_ratio~gene + codon  + aa + dataset_name, mean)
by_gene_plot <- merge(by_gene_plot, by_gene_plot_count, by=c("gene", "dataset_name"))
# make a label with sample numbers in it
by_gene_plot$labelled_name <- apply(by_gene_plot[,c("dataset_name", "genome_count")], 1, function(x)paste0(x[1], " (", x[2], ")"))
all_label <- unique(by_gene_plot$labelled_name)
all_label <- c(all_label[grep("bat", all_label)], all_label[grep("pangolin", all_label)], all_label[grep("charite", all_label)], all_label[grep("wuhan", all_label)])
by_gene_plot$labelled_name <- factor(by_gene_plot$labelled_name, levels=rev(all_label))

# reorder codon so they display in the order of blocks of codon coding for the same amino acid
by_gene_plot <- by_gene_plot[order(by_gene_plot$aa, by_gene_plot$codon),]
by_gene_plot$codon <- factor(by_gene_plot$codon, levels = unique(by_gene_plot$codon))
by_gene_plot$dataset_name <- factor(by_gene_plot$dataset_name, levels=rev(c(  "bat",  "pangolin",  "wuhan", "charite") ) )
x_axis_codon_order <- factor(unique(by_gene_plot$codon), levels=unique(by_gene_plot$codon))
x_axis_codon_colour <- factor(by_gene_plot$aa[!(duplicated(by_gene_plot$codon))], levels = unique(by_gene_plot$aa[!(duplicated(by_gene_plot$codon))]) )
levels(x_axis_codon_colour) <- rep(c("#fc8d62", "#377eb8"), length(levels(x_axis_codon_colour)))
x_axis_codon_colour <- as.character(x_axis_codon_colour)
# find the joints for change in aa
vline_positions = vector()
current_col=x_axis_codon_colour[1]
for(idx in 1: length(x_axis_codon_colour)){
	if(current_col != x_axis_codon_colour[idx]){
		vline_positions <- c(vline_positions, idx-0.5)
		current_col <- x_axis_codon_colour[idx]
	}
}

#mybreaks <- c(seq(0,2.4,0.3))
mybreaks <- c(seq(0,0.8,0.2))
include_genes <- c("multiple genes*", "ORF1a", "S", "N", "ORF3a", "M", "ORF8", "ORF7a", "E", "ORF6", "ORF7b", "ORF10")
by_gene_plot <- by_gene_plot[by_gene_plot$gene %in% include_genes, ]
by_gene_plot$gene <- factor(by_gene_plot$gene, levels = include_genes)
ggplot(by_gene_plot, aes(x=codon, y=labelled_name, fill=codon_ratio)) + geom_tile() + facet_wrap(~gene, ncol=1,scales = "free_y", strip.position='right') + 
	  theme_classic() +
	  theme(
		axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = x_axis_codon_colour),
		#text = element_text(size=6),
		strip.placement = "outside",
		strip.text.y = element_text(angle = 0),
		legend.direction = "horizontal",
		legend.position = "top",
		legend.key.width = unit(x = 2, units = "cm") ) +
	scale_fill_fermenter(palette = "Spectral",direction = -1,name="synonymous\ncodon ratio",limits=c(0,2.5), na.value = "grey50", breaks = mybreaks ) +
	ylab("") + xlab("") + geom_vline(xintercept=vline_positions, col='white', size=1)

ggsave("output/codon_usage_ratio.png", width=300, height=200,units='mm')



if(FALSE){

	## plot codon usage by amino acid
	#repeat grey colours
	full_df_list <- split(full_df, full_df$Ensembl_Gene)
	lapply(names(full_df_list), FUN=function(x)plot_codon_bias(full_df_list[[x]],x))
	plt <- plot_codon_bias(full_df, genome_collpase)
	aa_table_per_codon <- aa_table_per_codon[order(aa_table_per_codon$codon),]
	tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),  
						core=list(fg_params=list(fontsize=10)))
	tbl <- tableGrob(aa_table_per_codon, rows=NULL, theme=tt)
	g1 <- tableGrob(aa_table_per_codon[1:22,], rows=NULL, theme=tt)
	g2 <- tableGrob(aa_table_per_codon[23:43,] , rows=NULL, theme=tt)
	g3 <- tableGrob(aa_table_per_codon[44:64,] , rows=NULL, theme=tt)

	hcombined <-gtable_combine(g1,g2,g3)
	# Plot chart and table into one object
	png(paste(output_dir, "/multi_genes/", genome_collpase, ".png", sep=""), width=300,height=80, units='mm', res=300)
	grid.arrange(plt, hcombined,
				 ncol=2,
				 as.table=TRUE)
	dev.off()

}

