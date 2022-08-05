# proTSScall
Call active, dominant transcription start sites from PRO-seq data

	usage: Rscript proTSScall.R <tss.tes_list> <read_threshold>

	tss.tes_list	full path to TSS/TES list (e.g. hg38.basic_formakeheatmap.txt)
	read_threshold	TSS-proximal read cutoff, above which TSSs are considered active

Execute this script from the top level directory of NTC pipeline output (i.e. outdir). Running it will accomplish the following:

1. Organize all 1nt, 3' bedGraphs in bedGraphs/ into forward and reverse directories, in which all {forward,reverse} tracks are merged into a single {forward,reverse} bedGraph
2. makeheatmap is run on the merged bedGraphs, from TSS +/- 2kb, binning at 25bp
3. TSS-proximal read counts (TSS to +150) are calculated, and "inactive" TSSs are determined according to 'read_threshold'
4. From the remaining "active" TSSs, 1 dominant TSS is determined per gene (via TSS to +150 read counts), with ties won by the more upstream (strand-aware) TSS
5. Dominant TSSs are deduplicated (duplicates have same TSS start position) by keeping the longer transcript, with ties won by the lower ENSG number

Outputs (in pro_tss/):

	*_inactive_formakeheatmap.txt		inactive TSSs
	*_non-dominant_formakeheatmap.txt	active, non-dominant TSSs
	*_dominant_formakeheatmap.txt		dominant TSSs
	*_+250.TES.min400_formakeheatmap.txt	dominant gene bodies (min transcript length filter of 400)

---------------------------------	
Sample bash script:

	#!/bin/bash
	#SBATCH -p short
	#SBATCH -t 0-3
	#SBATCH --mem=50G
	#SBATCH -c 6
	
	module load gcc/6.2.0 R/4.0.1
	
	Rscript /path/to/proTSScall.R /path/to/tss_list 4
