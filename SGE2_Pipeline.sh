
#This pipeline is an example used to show the order in which the scripts are executed to analyze the data.

	#primers: X22 is F(1), R2 pair and X18 is the F3 R3 pair in this file on 7/21 (50%)
seq_type="NS" #or 
seq_dir="/net/shendure/vol9/seq/NEXTSEQ/170721_NS500488_0425_AHKVC5AFXX"
out_dir="/net/shendure/vol10/projects/SGE/nobackup/BRCA1/20170722_BRCA1_SGE2_X4X5X16X17X19X21X23"
sample_sheet="/net/gs/vol1/home/gf2/Documents/20170721_NS_sample_sheet_BRCA1_SGE2_X4X5X16X17X18X19X20X21X22X23.csv"
#splits first on comma, then +
cigar_comparisons="X4rL41_pre+X4rL41_post,X4rL42_pre+X4rL42_post,X5rL41_pre+X5rL41_post,X5rL42_pre+X5rL42_post,X19rL41_pre+X19rL41_post,X19rL42_pre+X19rL42_post,X21rL41_pre+X21rL41_post,X21rL42_pre+X21rL42_post,X23rL41_pre+X23rL41_post,X23rL42_pre+X23rL42_post,X16rL41_pre+X16rL41_post,X16rL42_pre+X16rL42_post"
amplicon_list="X4+X5+X16+X19+X21+X23"
##first gets split on +'s for each experiment, and then split on commas for key and pairings.
exp_groupings="X4rL41,X4rL41_pre,X4rL41_post,X4_lib,X4_neg,X4rL41_rnagDNA+X4rL42,X4rL42_pre,X4rL42_post,X4_lib,X4_neg,X4rL42_rnagDNA+X5rL41,X5rL41_pre,X5rL41_post,X5_lib,X5_neg,X5rL41_rnagDNA+X5rL42,X5rL42_pre,X5rL42_post,X5_lib,X5_neg,X5rL42_rnagDNA+X19rL41,X19rL41_pre,X19rL41_post,X19_lib,X19_neg,X19rL41_rnagDNA+X19rL42,X19rL42_pre,X19rL42_post,X19_lib,X19_neg,X19rL42_rnagDNA+X21rL41,X21rL41_pre,X21rL41_post,X21_lib,X21_neg,X21rL41_rnagDNA+X21rL42,X21rL42_pre,X21rL42_post,X21_lib,X21_neg,X21rL42_rnagDNA+X23rL41,X23rL41_pre,X23rL41_post,X23_lib,X23_neg,X23rL41_rnagDNA+X23rL42,X23rL42_pre,X23rL42_post,X23_lib,X23_neg,X23rL42_rnagDNA+X16rL41,X16rL41_pre,X16rL41_post,X16_lib,X16_neg,X16rL41_rnagDNA+X16rL42,X16rL42_pre,X16rL42_post,X16_lib,X16_neg,X16rL42_rnagDNA"
#### end variables section
echo "Variables defined."

mkdir $out_dir
mkdir $out_dir/fastq
module load bcl2fastq/2.16
bcl2fastq -R $seq_dir -o $out_dir/fastq --no-lane-splitting --sample-sheet $sample_sheet
cd $out_dir/fastq #this only changes the reference in which the rest of the script is interpreted...
python ~/bin/run_rename_no_Ss.py 
sh run_rename_no_Ss.sh
zgrep -c @$seq_type *R1_* | tee read_counts.txt
rm Undetermined*
module load fastqc/0.11.3
mkdir fastqc_out
fastqc *.fastq.gz -o fastqc_out &
echo "Running fastqc."
mkdir Seqprep
mkdir Seqprep/R1
mkdir Seqprep/R2
mkdir Seqprep/merged
python ~/bin/run_seqprep_BRCA1_pipeline.py
#editing this to wait longer
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
echo "Seqprep done and fastqc done."
cd Seqprep/merged
zgrep -c @$seq_type *.fastq.gz >> seqprep_read_counts.txt
mkdir no_Ns
python ~gf2/bin/run_remove_n_bases.py #operates with getcwd()
sh run_remove_n_bases.sh
echo "remove_n_bases done."
cd no_Ns
python ~/bin/run_cDNA_to_gDNA_SGE_pipeline.py
sh run_cDNA_to_gDNA.sh 
#a script to work on all samples with RNA in the name, to check for perfect ends flanking exons and if they are present, replace with intronic sequence using a hard-coded reference file specifying each amplicon's 5'/3' cDNA ends, full length, and 5'/3' intronic flanks
mkdir sam
python ~/bin/run_needle_to_sam_BRCA1_pipeline.py #switches from dashes to underscores, and uses a hardcoded reference file location with all 13 BRCA1 amplicons.
echo "Running needleall."
sh run_needle_to_sam.sh
wait
echo "Finished running needleall."
cd sam

mkdir cigar_counts
python ~/bin/20170330_BRCA1_SGE_pipeline_cigar_analyzer.py $cigar_comparisons #not informative for RNA based on full length processing
head -15 cigar_counts/*.txt >> cigar_counts/combined_cigar_counts.txt
mkdir variant_counts_no_thresh
#must rerun from here
python ~/bin/sam_to_edits_wRNA_pipeline_20170517.py $amplicon_list $exp_groupings
#writes minimally thresholded output (set in hardcode (0.000001 in any) to variant_counts folder
cd variant_counts_no_thresh
#ALL CADD SCORES NOW IN SAME DIRECTORY FOR BRCA1 SGE project and scripts point there.
mkdir final
python ~/bin/annotate_variants_BRCA1_pipeline_20170622.py
cd final
head -6 *summary.txt >> combined_editing_data.txt
echo 'Pipeline finished, variants aligned, counted, and annotated.'
#to run on old directory after running this script:

sh 20170629_SGE2_Pipeline_X2X3X15X18X20X22_patch_from_sam_with_new_adds.sh