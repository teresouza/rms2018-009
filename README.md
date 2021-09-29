# rms2018-009
This repository contains code used to process and generate plots for the rhabdomyosarcoma paper.
Two main folder (rnaSeq and wgs) contain the scripts that were used to process, filter and generate plots from rna-sequencing and whole genome sequencing, respectively.

## WGS folder content
* ```wgs_proc_RMS_call.R```: has paths and main function calls to process outputs from the GATK somatic copy number (called copy number alterations) and somatic single nucleotide variant (.vcf files). It sources ```wgs_proc_RMS_call.R``` that contains all user-defined functions
* ```wgs_proc_RMS_call.R```: contains user-defined functions to process the datasets. It sources ```wgs_circos_plot.R```, which generates circos plots for germline and tumor samples.
* ```wgs_freq_plots```, a standalone script that generates copy number frequency plots. 
* ```wgs_SNV_plots.R```, a standalone script that generates SNV plots.

## RNASeq folder content
* ```rnaSeq_proc.R```: contains the code used to process the RNA-sequencing count tables.
* ```rnaSeq_plots.R```: contains the code used to generate RNA-sequencing plots (correlogram heatmap and PCAs)
