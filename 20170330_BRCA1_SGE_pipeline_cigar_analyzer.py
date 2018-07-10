#20170314_BRCA1_cigar_analyzer.py

#Hardcoded file paths

#Pseudo #make dictionaries to count unique cigar strings w/i each file
        #compare top cigar strings from one file to top cigar strings from the other file 
        #spit out counts to simple txt files



import os
import sys
import subprocess

#function to generate a lookup table of cigar counts from a standard sam file w/ a 2 line header and cigars in index 5.
def build_cigar_dict(sam_file):
    line_count = 0
    cigar_counts = {}
    for line in sam_file:
        if line_count < 2:
            line_count += 1
            continue
        else:
            line_count += 1
            #pull out cigar string
            cigar_string = line.split('\t')[5]
            #print cigar_string
            if cigar_string in cigar_counts:
                cigar_counts[cigar_string]+=1
            else:
                cigar_counts[cigar_string]=1
    return cigar_counts #returns the final dictionary after parsing the file

def sort_counts_dict(counts_dict):
    return sorted(counts_dict.keys(), key=lambda k: -1*int(counts_dict[k]))

def get_read_count(counts_dict):
    total_reads = 0
    for my_key in counts_dict:
        total_reads += counts_dict[my_key]
    return total_reads

def write_top_100_outfile(dict1,dict2,file_path):
    out_file = open(file_path,'w')
    out_file.write('cigar\ts1_reads\ts2_reads\ts1_rpt\ts2_rpt\ts2_s1_ratio\n')
    outlines = 0
    sorted_dict1_keys = sort_counts_dict(dict1)
    sorted_dict2_keys = sort_counts_dict(dict2)
    total_dict1_reads = get_read_count(dict1)
    total_dict2_reads = get_read_count(dict2)
    sorted_dict1_key_count = len(sorted_dict1_keys)
    max_lines = 100
    if sorted_dict1_key_count < 100:
        max_lines = sorted_dict1_key_count
    while outlines < max_lines:
        outlines += 1
        #calculate all those things above, make them strings, and write them out after tab-joining.
        my_key = sorted_dict1_keys[outlines-1]
        d5_reads = dict1[my_key]
        if my_key in dict2:
            d11_reads = dict2[my_key]
        else:
            d11_reads = 0
        d5_rpt = float(d5_reads)/total_dict1_reads*1000
        d11_rpt = float(d11_reads)/total_dict2_reads*1000
        if d5_rpt != 0:
            d11_d5_ratio = d11_rpt/d5_rpt
        else:
            d11_d5_ratio = 9999
        output_list = [str(my_key),str(d5_reads),str(d11_reads),str(d5_rpt),str(d11_rpt),str(d11_d5_ratio)]
        out_file.write('\t'.join(output_list)+'\n')
    out_file.close()

#list of particular commands to analyze this data:
working_dir = os.getcwd()

#put all information into a large dictionary of dictionary, where each sam file is a key that points to it's counts dictionary
counts_dod = {}
for i in os.listdir(working_dir):
    if i.endswith('.sam'):
        print i
        counts_dod[str(i)] = build_cigar_dict(open(working_dir+'/'+i,'r'))

#uses a list of sample pairings for comparison as input sys.argv[1]
comp_list = sys.argv[1] #format the comparison list like: X15_lib+X15_lib,X17_lib+X17_lib...etc
comps = comp_list.split(',')
print comps
for comp in comps:
    samples = comp.split('+')
    s1 = samples[0]
    print s1
    s2 = samples[1]
    print s2
    print 'Comparing cigar counts for', s1, 'and', s2
    write_top_100_outfile(counts_dod[s1+'.sam'],counts_dod[s2+'.sam'],working_dir+'/cigar_counts/'+s1+'_'+s2+'_top100.txt')







