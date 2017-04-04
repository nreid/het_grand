# het_grand

many analyses are run using "parallel" on the UConn HPC cluster. 

directory 'align':
	* grandis_align.sh - base script that aligns and sorts individual Fundulus grandis samples
	* het_align.sh - base script that aligns and sorts individual Fundulus heteroclitus samples
	* $POP_align.slurm - these files invoke *_align.sh using 'parallel' for each population. submitted to slurm job scheduler
	* other scripts compress and sort split and discordant read files

directory 'genome_stats':
	* basecomp.slurm - calculates base composition in windows for reference genome
	* $POP_cov.slurm - calculates coverage in windows for each population or species
	* grandis_subpops_cov.slurm - loops through all grandis populations to calculate coverage

directory 'variants':
	* call_variants_subregion.sh - basic variant calling script. invoked using parallel
	* directory 'chunks':
		-300 scripts which each invoke call_variants_subregion.sh on a set of regions adding up to approximately 3.4mb