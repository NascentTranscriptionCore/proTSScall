# proTSScall
Call active, dominant transcription start sites from PRO-seq data

Usage: Rscript proTSScall.R <tss_list> <read_threshold>

	tss_list	full path to TSS list in makeheatmap format
	read_threshold	TSS-proximal read cutoff, *above* which TSSs are considered significantly active

Execute this script from the top level directory of NTC pipeline output (i.e. the parent directory of bedGraphs/, matrix/, etc). Running it will accomplish the following:

	1. Organize all *_{F,R}.bedGraph in bedGraphs/ into forward and reverse directories, in which all {F,R} samples are merged into a single {F,R} bedGraph
	2. makeheatmap is run on the merged bedGraphs, from TSS +/- 2kb, binning at 25bp (output in matrix/) 
	3. TSS-proximal read counts (TSS to +150) are calculated, and "inactive" TSSs are determined according to 'read_threshold'
	4. From the remaining "active" TSSs, 1 dominant TSS is determined per gene (via TSS to +150 read counts), with ties won by the more upstream TSS (strand-aware)
	5. Dominant TSSs are deduplicated (duplicates have same TSS start position) by keeping the longer transcript, with ties won by the lower ENSG number

Outputs (in pro_tss/):

	*_inactive_formakeheatmap.txt	inactive TSSs
	*_non-dominant_formakeheatmap.txt	active, yet non-dominant, TSSs
	*_dominant_formakeheatmap.txt	deduplicated, dominant TSSs
	*_tss.saf		TSS-proximal intervals (TSS to +150) for deduplicated dominant TSSs, in SAF format
	*_whole_gene.saf	whole gene intervals (TSS to transcript end) for deduplicated dominant TSSs, in SAF format
	*_gene_body.saf		gene body intervals (TSS+250 to TSS+2250) for deduplicated dominant TSSs, in SAF format

---------------------------------	
Sample bash script:

	#!/bin/bash
	#SBATCH -p short
	#SBATCH -t 0-3
	#SBATCH --mem=25G
	#SBATCH -c 6
	
	Rscript /path/to/proTSScall.R /path/to/tss_list 4

