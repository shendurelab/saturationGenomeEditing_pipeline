#run_cDNA_to_gDNA_SGE_pipeline.py
#updated on 6/30/2017 to work with 'rna' in file name and not 'RNA' 
import os
import sys
import subprocess
#hard-coded directories to point to python script, which has hardcoded directory for BRCA1_cDNA_data.txt
#will output a script to be run in shell

working_dir = os.getcwd()
shell_file = open('run_cDNA_to_gDNA.sh', 'w')
shell_file.write('cd '+working_dir+'\n')

for i in os.listdir(working_dir):
    if (i.endswith(".fastq")) and ('rna' in i): #will only work on fastq's that contain rna in sample name
        print i
        index_of_first_dot = i.find('.')
        index_of_extension = i.find('.fastq')
        index_of_first_r = i.find('r')
        amp = i[:index_of_first_r]
        sample_name = i[:index_of_first_dot] #followed by .merged.fastq
        before_extension_name = i[:index_of_extension]
        shell_file.write("python /net/gs/vol1/home/gf2/bin/cDNA_to_gDNA_SGE_pipeline.py "+i+' '+sample_name+'gDNA.merged.fastq ' +amp+' &\n')
    else:
        pass

shell_file.write('wait\n')

for i in os.listdir(working_dir):
    if (i.endswith(".fastq")) and ('rna' in i) and ('gDNA' not in i):
        shell_file.write('rm '+i+'\n')

shell_file.close()


#script format:  python