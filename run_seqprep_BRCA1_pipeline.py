#run_seqprep_BRCA1_pipeline.py
#to be run in the folder for the experiment with all the reads split by sample as fastq.gz files
#generates a shell script to run
#requires no underscores to be in sample names (later script will switch from dashes to underscores)
#outputs the stats from each seqprep call to a file called seqprep_stats.txt

import os
import sys
import subprocess

working_dir = os.getcwd()
print working_dir
command_file_name = working_dir+"/run_seqprep.sh"
command_file = open(command_file_name, 'w')

for i in os.listdir(os.getcwd()):
    if i.endswith("R1_001.fastq.gz"): 
        index_of_first_under = i.find('_')
        index_of_R1 = i.find('R1')
        sample_name = i[:index_of_first_under]
        file_name = i[:index_of_R1]
        # here are adapters for F and R seq primers
        command_file.write(r'/net/gs/vol1/home/gf2/bin/SeqPrep/./SeqPrep -f '+file_name+r'R1_001.fastq.gz -r '+file_name+r'R2_001.fastq.gz -1 Seqprep/R1/'+sample_name+r'.R1.fastq.gz -2 Seqprep/R2/'+sample_name+r'.R2.fastq.gz -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG -M 0.1 -s Seqprep/merged/'+sample_name+'.merged.fastq.gz -m 0.001 -q 20 -o 20 &\n')
    else:
        pass
command_file.write("wait\n")
command_file.close()


#notes for seqprep
'''
adapter trimming for A and B param's:
-A GGTTTGGAGCGAGATTGATAAAGT #if PU1R is used (as seen in R1):  
-B CTGAGCTCTCTCACAGCCATTTAG #if PU1L is used (as seen in R2): 
-M 0.1
'''