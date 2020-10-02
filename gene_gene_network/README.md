# Coronahack2020 gene-gene network analysis

This pipeline takes the Prokka output (from the host-data folder) and Ensembl v100 ASM985889v3 (in the input_ref_genes folder) and generates gene-gene network analysis as .gml files.

## Usage:
> ./generate_network.sh

This would create all_nuc.gml (network graph using DNA sequence) and all_aa.gml (network graph using amino acid sequence). These files can be opened using Graphia (https://graphia.app/). 

## Requirements:
- blast+ (version 2.7.1)
- python (version 3.4)
