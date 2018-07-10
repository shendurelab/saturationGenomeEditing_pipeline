README_BRCA1_SGE_custom_scripts.txt

This README file describes BRCA1 SGE data processing scripts.

All sequencing was performed on an Illumina NextSeq or MiSeq machine, and processed using bcl2fastq v2.16 from Illumina. The resulting fastq.gz files are available for download.

Code is provided for documentation of the analysis that was performed. Please see the Methods section of the paper for additional details.

Programs used for analysis:
bcl2fastq v2.16
Python  2.7.3
SeqPrep (available at https://github.com/jstjohn/SeqPrep)
fastqc v0.11.3
EMBOSS v6.4.0
R v3.1.3
RStudio v1.0.153

All scripts are located in the folder here:  ~/SGE/BRCA1/scripts_for_sharing

Custom python scripts used in analysis are as follows:

run_seqprep_BRCA1_pipeline.py
run_remove_n_bases.py
run_cDNA_to_gDNA_SGE_pipeline.py
run_needle_to_sam_BRCA1_pipeline.py
20170330_BRCA1_SGE_pipeline_cigar_analyzer.py
sam_to_edits_wRNA_pipeline_20170517.py
annotate_variants_BRCA1_pipeline_20170622.py

A pipeline used to call these scripts in the correct order is here:

SGE2_Pipeline.sh

Annotated, tab-delimited text files were imported to R Markdown for further analysis.

R_markdown scripts used to perform the analysis are as follows:

BRCA1_SGE1_X_all_ggplot_removed_20171012.Rmd
positional_modeling_SGE1_v6_20171013_noggplot.Rmd
BRCA1_SGE1_global_20171119_noggplot.Rmd

BRCA1_SGE2_X_all_ggplot_removed_20170801.Rmd
positional_modeling_SGE2_v6_20171012_noggplot.Rmd
BRCA1_SGE2_global_20180111.Rmd 

BRCA1_SGE_SxRepro_correlations_20171215.Rmd
BRCA1_SGE1_indel_analysis20171218.Rmd












