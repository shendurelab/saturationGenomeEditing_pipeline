
#annotate_variants_BRCA1_pipeline_20170622.py
#also imports cadd scores with exac data.
#current patch to fix variant annotation in the codon's of mutated PAMs
#updated version takes RNA sample as 5th in exp_grouping
#has more CADD annotation columns including transcript information
#specifies the CCDS from which to draw CADD information
#will provide 'REF' annotation for all CADD entries in reference alleles instead of 'NA'
#calls error rates based on the neg sample.
#gives the pSNV-ref alleles a 'pos' to guarantee that they get plotted in seq-func maps (make REF a consequence for coloring by cons. - black)
#summary output includes how many programmed mutations (i.e. oligos that were ordered) are in each library, and how many possible
#ClinVar annotations from clinvar 1/02/2018 -- located here:
    #/net/shendure/vol10/projects/SGE/nobackup/BRCA1/ClinVar_BRCA1_values.txt
    #updated from previous 5/22/17 version used up until 1/2/18.
    #One duplicate annotation removed manually from excel file. Had to copy out of excel into sublime text and save as .txt to get proper tdl format for reading into python (new line character issues)
#also includes 'nudge' values to adjust position on plots to show 3 variants per site.

import sys
import os
import subprocess
import time
from datetime import datetime
startTime = datetime.now()
#CRITICAL -- this line adds to the path the location of your binary folder so you can import modules from there (in this case, dict_tools_11192015.py)
sys.path.insert(0,'/net/gs/vol1/home/gf2/bin')
import dict_tools_11192015 as dict_tools


#--------------------------------------------------------------------------------------------------------------------------------------------
#HERE ARE BRCA1 SGE SPECIFIC CONDITIONS HARDCODED:
my_dir =  os.getcwd()
output_subdir = '/final'
working_dir = os.getcwd()
if my_dir == working_dir:
    print 'annotate_variants script is running in', my_dir

#Common thing in files to be analyzed here:  
common_file_id = '.txt'
amplicon_list = "X2+X3+X4+X5+X15+X16+X17+X18+X19+X20+X21+X22+X23".split('+') #will pull in amplicon info for all of the exons.

#REFERENCE FILES (in fasta format -- takes second line of file), hardcoded to look in main BRCA1 fasta storage:
fasta_dir = '/net/shendure/vol10/projects/SGE/nobackup/BRCA1/fasta'
ref_seqs = {}
for amplicon in amplicon_list:
    ref_file = open(fasta_dir+'/'+amplicon+'.fa', 'r')
    ref_header = ref_file.readline()
    ref_seq = ref_file.readline().strip()
    ref_seqs[amplicon] = ref_seq
    ref_file.close()

#editing info for each experiment, including region that is mutated
#This can stay hardcoded for all BRCA1 experiments, now.
edits_per_amp = {}
editing_info = open('/net/shendure/vol10/projects/SGE/nobackup/BRCA1/BRCA1_editing_data.txt', 'r')
editing_info_header = editing_info.readline()
for line in editing_info:
    edits = line.strip().split()
    #each amplicon in the file is a key that points to a list containing the info for that amplicon [5' HDR, 3' HDR, cut_pos, cut_end, mut_start, mut_end] 
    edits_per_amp[edits[0]] = edits[1:]

#cadd file location -- note this has been updated to point to the exac-enhanced files w. the e.cadd :
cadd_file_dict = {}
for amplicon in amplicon_list:
    cadd_file_dict[amplicon] = open('/net/shendure/vol10/projects/SGE/nobackup/BRCA1/cadd/'+amplicon+'e.cadd','r') 

clinvar_dict = {}
#make a clinvar dict, and add a simpler pathogenicity annotation system to each cadd_lookup --> Clinical_Sign. call.
cv_line_count = 0
unique_clin_sigs = []
header_indices = {} #dictionary pointing from header title to it's index
my_clinvar_file = open('/net/shendure/vol10/projects/SGE/nobackup/BRCA1/ClinVar_BRCA1_values.txt', 'r')
for line in my_clinvar_file:
    line_values = line.strip().split('\t') #txt tab-delim file of scores
    if cv_line_count == 0:
        cv_line_count += 1
        for header in line_values:
            my_header = header.strip()
            header_indices[my_header] = line_values.index(header)
        print header_indices
    elif cv_line_count != 0:
        cadd_key = line_values[header_indices['CADD_Key']]
        clin_sig = line_values[header_indices['Clinical_Significance']]
        if clin_sig not in unique_clin_sigs:
            unique_clin_sigs.append(clin_sig)
        if clin_sig not in ['Benign','Likely benign','Uncertain significance','Likely pathogenic','Pathogenic','Conflicting interpretations of pathogenicity']:
            if clin_sig == '"Conflicting interpretations of pathogenicity, not provided"':
                simp_clin_sig = 'Conflicting interpretations of pathogenicity'
            elif clin_sig == '"Benign/Likely benign, not provided"':
                simp_clin_sig = 'Likely benign'
            elif clin_sig == '"Pathogenic/Likely pathogenic, not provided"':
                simp_clin_sig = 'Likely pathogenic'
            elif clin_sig == 'Benign/Likely benign':
                simp_clin_sig = 'Likely benign'
            elif clin_sig == 'Pathogenic/Likely pathogenic':
                simp_clin_sig = 'Likely pathogenic'
        else:
            simp_clin_sig = clin_sig
        if cadd_key in clinvar_dict:
            print 'error '+cadd_key+' duplicated in ClinVar.'
        clinvar_dict[cadd_key] = [clin_sig,simp_clin_sig]

#sequencing amplicon genomic location and orientation
seq_start_dict = {'X2':41276192,'X3':41267891,'X4':41258618,'X5':41257033,'X15':41223093,'X16':41219775,'X17':41216050,'X18':41215466,'X19':41209219,'X20':41203207,'X21':41201273,'X22':41199777,'X23':41197877,'BRCA1x18':41216050}
orientation = '-'
BRCA1_CCDS = 'CCDS11456.2'
#line_count is set up for a 2-line header (default for tabix pull of cadd file)

#--------------------------------------------------------------------------------------------------------------------------------------------
def get_snv_mismatch(edit_string,exp_hdr_edits): #will take an edit_string (with 1 mismatch only) and return the cadd_lookup_key
    edits = edit_string.split(',')
    mismatch_edits = []
    for edit in edits:
        if 'X' in edit:
            mismatch_edits+=[edit]
    copied_list = mismatch_edits[:] #have to copy the list so you're not iterating over a list that is changing based on removal.
    for mm_edit in copied_list:
        if mm_edit in exp_hdr_edits:
            mismatch_edits.remove(mm_edit)
        else:
            pass
    if len(mismatch_edits) == 1:
        snv_edit = mismatch_edits[0]
        return(snv_edit)
    else:
        print "Error: The edit edit_string", edit_string, 'did not contain a SNV.'

#--------------------------------------------------------------------------------------------------------------------------------------------
def mismatch_to_cadd_lookup(mismatch_data,chrom_index,orientation): #will take an edit_string (with 1 mismatch only) and return the cadd_lookup_key
    edit_data = mismatch_data.strip().split('-') #[0]=position
    pos = int(edit_data[0])
    #such that first base in amplicon == chrom_index
    base = ''
    if '-' in orientation:
        chr_pos = int(chrom_index)-(pos-1)
        base = edit_data[2]
        if base == 'A':
            base = 'T'
        elif base == 'T':
            base = 'A'
        elif base == 'C':
            base = 'G'
        elif base == 'G':
            base = 'C'
    elif '+' in orientation:
        chr_pos = int(chrom_index)+(pos-1)
        base = edit_data[2]
    return(str(chr_pos)+base)
#--------------------------------------------------------------------------------------------------------------------------------------------

def get_cadd_data_for_variant(cadd_lookup_key, cadd_data): #cadd_data includes all annotations, this will return only the ones you want:
    variant_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],cadd_data[cadd_lookup_key][10],cadd_data[cadd_lookup_key][12],cadd_data[cadd_lookup_key][110],cadd_data[cadd_lookup_key][111],cadd_data[cadd_lookup_key][114],cadd_data[cadd_lookup_key][115],cadd_data[cadd_lookup_key][18],cadd_data[cadd_lookup_key][22],cadd_data[cadd_lookup_key][24],cadd_data[cadd_lookup_key][25],cadd_data[cadd_lookup_key][29],cadd_data[cadd_lookup_key][38],cadd_data[cadd_lookup_key][81],cadd_data[cadd_lookup_key][82],cadd_data[cadd_lookup_key][85],cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],cadd_data[cadd_lookup_key][108],cadd_data[cadd_lookup_key][109],cadd_data[cadd_lookup_key][112],cadd_data[cadd_lookup_key][113],cadd_data[cadd_lookup_key][116],cadd_data[cadd_lookup_key][117],cadd_data[cadd_lookup_key][118],cadd_data[cadd_lookup_key][120],cadd_data[cadd_lookup_key][121]]

    return variant_CADD_data

#--------------------------------------------------------------------------------------------------------------------------------------------

def get_CDSpos_for_variant(cadd_lookup_key, cadd_data): #cadd_data includes all annotations, this will return only the 10 you want:
    CDSpos = cadd_data[cadd_lookup_key][98]
    if CDSpos.isdigit():
        CDSpos = int(CDSpos)
    return CDSpos

#--------------------------------------------------------------------------------------------------------------------------------------------

def header_index(header_line, category): #takes tab-delimited headers only, and category must be a string that defines a category
    header_list = header_line.strip().split('\t') #header list will have 'variant' category as well
    var_dict_index = header_list.index(category)-1
    return var_dict_index

#--------------------------------------------------------------------------------------------------------------------------------------------
def get_del_sizes(edit_string):
	#split the edit string on commas into a list
	del_sizes = []
	#This will turn the edits into a list of strings, regardless of how many there are (0 will be ['WT']...)
	edit_list = edit_string.split(',')
	#iterate over each list item and ask if there's a deletion
	for edit in edit_list:
		if '-D' in edit:
			index_of_D = edit.find('-D')
			del_sizes.append(int(edit[index_of_D+2:]))
	return del_sizes

#--------------------------------------------------------------------------------------------------------------------------------------------
def get_del_starts(edit_string):
	del_starts = []
	edit_list = edit_string.split(',')
	#iterate over each list item and ask if there's a deletion
	for edit in edit_list:
		if '-D' in edit:
			index_of_D = edit.find('-D')
			del_start = int(edit[:index_of_D])
			del_starts.append(del_start)
	return del_starts

#--------------------------------------------------------------------------------------------------------------------------------------------
def get_del_bases(edit_string):
	total_del_bases = []
	edit_list = edit_string.split(',')
	#iterate over each list item and ask if there's a deletion
	for edit in edit_list:
		if '-D' in edit:
			index_of_D = edit.find('-D')
			del_start = int(edit[:index_of_D])
			del_size = int(edit[index_of_D+2:])
			del_bases = range(del_start,del_start+del_size) 
			total_del_bases.append(del_bases)
	return total_del_bases #THIS IS A LIST of LISTS where each sub-list is all the bases deleted in a single deletion
	#DEL LOCATIONS ARE 1-indexed for the sequencing read -- neeed to be converted in sum_del_bases to 0-indexing

#--------------------------------------------------------------------------------------------------------------------------------------------
def sum_del_bases(running_totals,list_of_del_positions):
	#Example:  running_totals = [0,0,1,1,2,2,2,0,0,0], list_of...=[1,4,5,9]
	#Want to increment the list at indexes 0,3,4,8 to [1,0,1,2,3,2,2,0,1,0]...
	for del_po in list_of_del_positions:
		#first entry in running totals is the FIRST BASE
		running_totals[del_po-1]+=1
#-------------------------------------------------------------------------------------------------------------------------------------------

#FOR INSERTIONS
#-------------------------------------------------------------------------------------------------------------------------------------------
def get_insert_sizes(edit_string):
	insert_sizes = []
	edit_list = edit_string.split(',')
	for edit in edit_list:
		if '-I' in edit:
			index_of_I = edit.find('-I')
			index_of_second_dash = edit.rfind('-')
			insert_sizes.append(int(edit[index_of_I+2:index_of_second_dash]))
	return insert_sizes

#--------------------------------------------------------------------------------------------------------------------------------------------
def get_insert_starts(edit_string):
	ins_starts = []
	edit_list = edit_string.split(',')
	for edit in edit_list:
		if '-I' in edit:
			index_of_I = edit.find('-I')
			ins_start = int(edit[:index_of_I])
			ins_starts.append(ins_start)
	return ins_starts

#-------------------------------------------------------------------------------------------------------------------------------------------
def get_insert_bases(edit_string):
	total_ins_bases = []
	edit_list = edit_string.split(',')
	for edit in edit_list:
		if '-I' in edit:
			index_of_second_dash = edit.rfind('-')
			insert_bases = edit[index_of_second_dash+1:]
			total_ins_bases.append(insert_bases)
	return total_ins_bases #THIS IS A LIST of where each member is a single insertions sequence; if one insertion, then bases will be
#-------------------------------------------------------------------------------------------------------------------------------------------

def increment_variant_count (read_counts):  #will take a number and return 0 if absent and 1 if present
    if read_counts > 0:
        return 1
    else:
        return 0

def increment_variant_count_thresh (read_counts, threshold):  #will do same as above, but instead of being 1+ times, must be threshold or more.
    if read_counts >= int(threshold):
        return 1
    else:
        return 0

def mutagenize(my_string):
    upper_string = my_string.upper()
    #the first element will be the wild-type sequence
    all_variants = [upper_string]
    for i in range (0,len(upper_string)):
        for j in ("A", "C", "G", "T"):
            if j != upper_string[i]:
                all_variants += [upper_string[:i]+j+upper_string[i+1:]]
    return all_variants

def get_all_hdr_variants(amplicon,ref_seq_dict,edits_per_amp_dict): #requires mutagenize function above
    #here, you want to feed it amplicon, ref_seqs, and edits_per_amp as input
    ref_seq = ref_seq_dict[amplicon]
    all_hdr_variants = [ref_seq] #list will start with total wt
    my_editing_info = edits_per_amp_dict[amplicon]
    #lists of 5' and 3' edits
    edits5 = my_editing_info[0].strip().split(',')
    edits3 = my_editing_info[1].strip().split(',')
    cut_pos = int(my_editing_info[2])
    cut_end = my_editing_info[3]
    mut_start = int(my_editing_info[4])
    mut_end = int(my_editing_info[5])
    all_pam_edits = edits5 + edits3
    mut_ref_seq = str(ref_seq)
    for edit in all_pam_edits:
        edit_location = int(edit[:edit.find('-')])
        new_base = edit[-1]
        mut_ref_seq = mut_ref_seq[:edit_location-1]+new_base+mut_ref_seq[edit_location:]
    #mut_ref_seq now has all marker SNVs in it
    f_adapt = mut_ref_seq[:mut_start-1]
    mut_seq = mut_ref_seq[mut_start-1:mut_end]
    r_adapt = mut_ref_seq[mut_end:]
    all_mut_seqs = mutagenize(mut_seq)
    for seq in all_mut_seqs:
        all_hdr_variants.append(f_adapt+seq+r_adapt)
    return all_hdr_variants #complete list with wt wt, all cSNVs and 

def get_pos_nudge(ref_base,alt_base): #assumes gene is on the negative strand
    pos_nudge = 0
    CDSpos_nudge = 0
    rev_comp = ''
    if alt_base == 'A':
        rev_comp = 'T'
        pos_nudge = -0.25
        CDSpos_nudge = 0.25
    elif alt_base == 'T':
        rev_comp = 'A'
        pos_nudge = 0.25
        CDSpos_nudge = -0.25
    elif alt_base == 'C':
        rev_comp = 'G'
        if ref_base == 'A':
            pos_nudge = -0.25
            CDSpos_nudge = 0.25
        else:
            pos_nudge == 0
            CDSpos_nudge == 0
    elif alt_base == 'G':
        rev_comp = 'C'
        if ref_base == 'T':
            pos_nudge = 0.25
            CDSpos_nudge = -0.25
        else:
            pos_nudge = 0
            CDSpos_nudge = 0
    return [rev_comp,pos_nudge,CDSpos_nudge]

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def translate_codon(codon):
    return codon_table[codon]

def get_rev_comp(ref_base):
    ref_base = ref_base.upper()
    if ref_base == 'A':
        rc_base = 'T'
    elif ref_base == 'T':
        rc_base = 'A'
    elif ref_base == 'C':
        rc_base = 'G'
    elif ref_base == 'G':
        rc_base = 'C'
    return rc_base

#-------------------------------------------------------------------------------------------------------------------------------------------
#This can stay hardcoded for all BRCA1 experiments, now.
edits_per_amp = {}
editing_info = open('/net/shendure/vol10/projects/SGE/nobackup/BRCA1/BRCA1_editing_data.txt', 'r')
editing_info_header = editing_info.readline()
for line in editing_info:
    edits = line.strip().split()
    #each amplicon in the file is a key that points to a list containing the info for that amplicon [5' HDR, 3' HDR, cut_pos, cut_end, mut_start, mut_end] 
    edits_per_amp[edits[0]] = edits[1:]

#get the cadd scores out for each amplicon:
cadd_data_dod = {}
CDS_dict = {} #format this to have CDSpos (integer) point to the reference base (genome), coding base (transcript), and 'oAA' in list []
for amplicon in amplicon_list:
    cadd_data_dod[amplicon] = {}
    line_count = 0
    chr_start = "not found"
    cadd_file = cadd_file_dict[amplicon]
    for line in cadd_file:
        line_count += 1
        if line_count > 2:
            variant_data = line.strip().split('\t')
            if line_count == 3:
                #get the chr_start location from 
                chr_start = int(variant_data[1])
            variant_key = str(variant_data[1])+str(variant_data[4]) #chromosomal location + cadd_lookup data
            if variant_key not in cadd_data_dod[amplicon]: #adding first transcript annotation first.
                cadd_data_dod[amplicon][variant_key] = variant_data
                CCDS = variant_data[94]
                CDSpos = variant_data[98]
                if CCDS == BRCA1_CCDS:
                    if (CDSpos not in CDS_dict) and (CDSpos != 'NA'):
                        ref_base = variant_data[2]
                        transcript_base = ref_base
                        if orientation == '-': #will RC the ref base if transcript is on the minus strand.
                            transcript_base = get_rev_comp(ref_base)
                        CDS_dict[int(CDSpos)] = [ref_base,transcript_base,variant_data[107]]
            elif variant_key in cadd_data_dod[amplicon]:
                #check that the transcript annotation matches BRCA1, if it does, replace the key:
                CCDS = variant_data[94]
                CDSpos = variant_data[98]
                if CCDS == BRCA1_CCDS:
                    cadd_data_dod[amplicon][variant_key] = variant_data
                    if (CDSpos not in CDS_dict) and (CDSpos != 'NA'):
                        ref_base = variant_data[2]
                        transcript_base = ref_base
                        if orientation == '-': #will RC the ref base if transcript is on the minus strand.
                            transcript_base = get_rev_comp(ref_base)
                        CDS_dict[int(CDSpos)] = [ref_base,transcript_base,variant_data[107]]
                else:
                    pass
        else:
            pass

#summary dictionary:  format experiment --> summary data (will need a header, too)
summary_data = {}

for i in os.listdir(working_dir):
    if common_file_id in i:
        #retrieve the experiment name to create entry for summary stats file
        #exp
        exp = str(str(i)[:str(i).index('.')])
        #DEFINE THE AMPLICON HERE:
        if 'r' in exp:
            amplicon = exp[:exp.find('r')]
        else:
            amplicon = exp
        ref_seq = ref_seqs[amplicon]
        seq_start = seq_start_dict[amplicon]
        oligo_variants = get_all_hdr_variants(amplicon,ref_seqs,edits_per_amp)
        len_oligo_variants = len(oligo_variants)
        ####
        #Extract editing info for the experiment (hdr mm edits, and the mutated region's start and end)
        mut_end = edits_per_amp
        snv_edits_5 = edits_per_amp[amplicon][0]
        snv_edits_3 = edits_per_amp[amplicon][1]
        exp_hdr_edits = snv_edits_5.split(',')+snv_edits_3.split(',')
        hdr_codons = []
        cadd_data = cadd_data_dod[amplicon]
        exp_mut_start = int(edits_per_amp[amplicon][-2])
        exp_mut_end = int(edits_per_amp[amplicon][-1])
        pam_codon_data = {} #format this CDSpos(int) --> [index_in_codon)(0,1,2),hdr_o_codon,hdr_n_codon]
        all_pam_edit_locations = []
        edit_5_locations = []
        for my_5_edit in snv_edits_5.split(','):
            edit_5_locations.append(int(my_5_edit[:my_5_edit.find('-X')]))
        edit_3_locations = []
        for my_3_edit in snv_edits_3.split(','):
            edit_3_locations.append(int(my_3_edit[:my_3_edit.find('-X')]))



        for hdr_edit in exp_hdr_edits:
            hdr_edit_pos = int(hdr_edit.split('-')[0])
            all_pam_edit_locations.append(hdr_edit_pos)
            hdr_base = hdr_edit.split('-')[-1] #from amplicon reference, which is in frame.
            if (hdr_edit_pos < exp_mut_start-2) or (exp_mut_end+2 < hdr_edit_pos): #requires the hdr_edit position to be within 2 bp sge region.
                print hdr_edit_pos, 'hdr_snv out of mutated range from:', exp_mut_start, exp_mut_end
                continue
            else:
                #lookup the codon if it exists
                cadd_lookup_key = mismatch_to_cadd_lookup(hdr_edit,seq_start,orientation)
                hdr_edit_cadd_data = cadd_data_dod[amplicon][cadd_lookup_key]
                hdr_edit_CDSpos = hdr_edit_cadd_data[98]
                if hdr_edit_CDSpos.isdigit(): #only converts to an integer if the string is a + integer (digit)
                    hdr_edit_CDSpos = int(hdr_edit_CDSpos)
                    rb_tb_aa = CDS_dict[hdr_edit_CDSpos] #ref_base, transcript_base, aa
                    hdr_o_codon = '' #original codon
                    hdr_n_codon = ''
                    hdr_codon_pos = [] #list of CDSpositions defining codon
                    snv_pos_in_frame = hdr_edit_CDSpos%3
                    if snv_pos_in_frame == 0: #pos 0 is third base (evenly divisble CDS number...)
                        AApos = hdr_edit_CDSpos/3
                        snv_pos_in_frame = 3
                        hdr_codon_pos = [hdr_edit_CDSpos-2,hdr_edit_CDSpos-1,hdr_edit_CDSpos] #list of positions defining the codon
                        hdr_o_codon = CDS_dict[hdr_edit_CDSpos-2][1]+CDS_dict[hdr_edit_CDSpos-1][1]+CDS_dict[hdr_edit_CDSpos][1]
                        hdr_n_codon = CDS_dict[hdr_edit_CDSpos-2][1]+CDS_dict[hdr_edit_CDSpos-1][1]+hdr_base
                        print amplicon, translate_codon(hdr_o_codon), CDS_dict[hdr_edit_CDSpos][2], 'match? code correct...'
                        print hdr_base, hdr_o_codon, hdr_n_codon, translate_codon(hdr_o_codon), translate_codon(hdr_n_codon), 'match? pam snv is syn...'
                        for j in range(0,3):
                            pam_codon_data[hdr_codon_pos[j]] = [j,hdr_o_codon,hdr_n_codon]
                    elif snv_pos_in_frame == 1:
                        AApos = hdr_edit_CDSpos/3 + 1
                        hdr_o_codon = CDS_dict[hdr_edit_CDSpos][1]+CDS_dict[hdr_edit_CDSpos+1][1]+CDS_dict[hdr_edit_CDSpos+2][1]
                        hdr_n_codon = hdr_base+CDS_dict[hdr_edit_CDSpos+1][1]+CDS_dict[hdr_edit_CDSpos+2][1]
                        hdr_codon_pos = [hdr_edit_CDSpos,hdr_edit_CDSpos+1,hdr_edit_CDSpos+2] #list of CDS positions defining the codon
                        for j in range(0,3):
                            pam_codon_data[hdr_codon_pos[j]] = [j,hdr_o_codon,hdr_n_codon]
                    elif snv_pos_in_frame == 2:
                        AApos = hdr_edit_CDSpos/3 + 1
                        hdr_o_codon = CDS_dict[hdr_edit_CDSpos-1][1]+CDS_dict[hdr_edit_CDSpos][1]+CDS_dict[hdr_edit_CDSpos+1][1]
                        hdr_n_codon = CDS_dict[hdr_edit_CDSpos-1][1]+hdr_base+CDS_dict[hdr_edit_CDSpos+2][1]
                        hdr_codon_pos = [hdr_edit_CDSpos-1,hdr_edit_CDSpos,hdr_edit_CDSpos+1] #list of CDS positions defining the codon
                        for j in range(0,3):
                            pam_codon_data[hdr_codon_pos[j]] = [j,hdr_o_codon,hdr_n_codon]

                else: #what to do if the PAM SNV is non-coding? nothing for now -- create a flag to indicate SNVs within 2 bp of PAM edit later.
                    pass
        print pam_codon_data
        # at this point for each pam edit is in the mutated region, and is covered with SNVs, you have a lookup table to get the codons and the postion in the codons by CDSpos of each variant (below)
        #now adapt code below to see if CDSpos of each variant is in the lookup table, and if so, to change the entries for AA, and conseq.
        #get the data in from the file
        header_and_dict = dict_tools.import_matrix2dict(i)
        header_list = header_and_dict[0]
        var_dict = header_and_dict[1]
        #make a list of the order the files are in to preserve it.
        ordered_variants = []
        line_count = 0
        for line in open(i,'r'):
            #skip the header
            if line_count != 0:
                ordered_variants.append(line.strip().split('\t')[0]) #add variant to list
            else:
                line_count += 1
        #HERE ARE THE HEADERS FROM PREVIOUS SCRIPT (for indexing)
        header_string = 'variant\tcigar\text_cigar\tedit_string\tpre\tpost\tlib\tneg\trna\tpre_pseudo\tpost_pseudo\tlib_pseudo\tneg_pseudo\trna_pseudo\tpre_freq\tpost_freq\tlib_freq\tneg_freq\trna_freq\tpre_pseudo_freq\tpost_pseudo_freq\tlib_pseudo_freq\tneg_pseudo_freq\trna_pseudo_freq\tpre_lib_ratio\tpost_pre_ratio\tpost_lib_ratio\trna_pre_ratio\trna_post_ratio\t'
        header_string+='hdr5\thdr3\t'
        header_string+='mismatch_count\tdel_count\tins_count\t'
        header_string+='all_hdr_variant\t'
        header_string+='pHDR_pre\tpHDR_post\tpHDR_lib\tpHDR_neg\tpHDR_rna\ttHDR_pre\ttHDR_post\ttHDR_lib\ttHDR_neg\ttHDR_rna\t'
        header_string+='tHDR_pre_pseudo_freq\ttHDR_post_pseudo_freq\ttHDR_lib_pseudo_freq\ttHDR_neg_pseudo_freq\ttHDR_rna_pseudo_freq\ttHDR_pre_lib_ratio\ttHDR_post_pre_ratio\ttHDR_post_lib_ratio\ttHDR_rna_pre_ratio\ttHDR_rna_post_ratio\t'

        if len(var_dict[ordered_variants[0]]) == (len(header_string.strip().split('\t'))-1):
            print 'Length of header matches length of data columns.'
        else:
            print 'Lengths of header list and variant dictionary are not equal!'

        #indices of key variaint identifiers
        index_pre = header_index(header_string,'pre')
        index_post = header_index(header_string,'post')
        index_lib = header_index(header_string,'lib')
        index_neg = header_index(header_string,'neg')
        index_rna = header_index(header_string,'rna')
        index_hdr5 = header_index(header_string,'hdr5')
        index_hdr3 = header_index(header_string,'hdr3')
        index_mismatch_count = header_index(header_string,'mismatch_count')
        index_del_count = header_index(header_string,'del_count')
        index_ins_count = header_index(header_string,'ins_count')
        index_edit_string = header_index(header_string,'edit_string')
        index_all_hdr_variant = header_index(header_string,'all_hdr_variant')

        #Summary stats to keep track of, write to an experimental summary TDL file at end?
        ref_len = len(ref_seq)
        #for calculating total HDR rate, initialize a summary data link -- use raw read counts here (no need for ratios)
        summary_header = '\t \n'
        summary_data[exp] = {}
        for timepoint in ['pre','post','lib','neg','rna']:

            summary_data[exp][timepoint] = {}
            summary_data[exp][timepoint]['total_reads'] = 0
            summary_data[exp][timepoint]['c_hdr_variants'] = 0
            summary_data[exp][timepoint]['c_hdr_reads'] = 0
            summary_data[exp][timepoint]['c_hdr_wt_variants'] = 0
            summary_data[exp][timepoint]['c_hdr_wt_reads'] = 0
            summary_data[exp][timepoint]['c_hdr_snv_variants'] = 0
            summary_data[exp][timepoint]['c_hdr_snv_reads'] = 0
            summary_data[exp][timepoint]['p_hdr_variants'] = 0
            summary_data[exp][timepoint]['p_hdr_reads'] = 0
            summary_data[exp][timepoint]['t_hdr_variants'] = 0
            summary_data[exp][timepoint]['t_hdr_reads'] = 0
            summary_data[exp][timepoint]['p_hdr_5_snv_reads'] = 0
            summary_data[exp][timepoint]['p_hdr_5_snv_variants'] = 0
            summary_data[exp][timepoint]['p_hdr_5_wt_reads'] = 0
            summary_data[exp][timepoint]['p_hdr_5_wt_variants'] = 0
            summary_data[exp][timepoint]['p_hdr_3_snv_reads'] = 0
            summary_data[exp][timepoint]['p_hdr_3_snv_variants'] = 0
            summary_data[exp][timepoint]['p_hdr_3_wt_reads'] = 0
            summary_data[exp][timepoint]['p_hdr_3_wt_variants'] = 0
            summary_data[exp][timepoint]['full_length_variants'] = 0
            summary_data[exp][timepoint]['full_length_reads'] = 0
            summary_data[exp][timepoint]['wt_wt_variants'] = 0
            summary_data[exp][timepoint]['wt_wt_variants_reads'] = 0
            summary_data[exp][timepoint]['wt_snv_variants'] = 0
            summary_data[exp][timepoint]['wt_snv_reads'] = 0
            summary_data[exp][timepoint]['oligo_variants'] = 0
        #Keep track of variants in key categories to normalize:
        c_hdr_snv_variants = []
        c_hdr_wt_variants = []
        
        wt_snv_variants = []
        wt_wt_variants = []


        add_to_header_1 = 'c_hdr\tp_hdr_5\tp_hdr_3\tno_hdr\tno_indels\tsingle_snv\tno_snv\twt_wt\twt_snv\tc_hdr_wt\tc_hdr_snv\tp_hdr_5_wt\tp_hdr_5_snv\tp_hdr_3_wt\tp_hdr_3_snv\twt_single_del\tc_hdr_single_del\t'
        add_to_header_2 = '\t'.join(['chr', 'pos', 'ref', 'alt', 'conseq', 'consdetail', 'polyphen-cat', 'polyphen-val', 'CADD-raw', 'CADD-phred','priPhCons', 'mamPhyloP', 'GerpN', 'GerpS', 'mutIndex', 'fitCons', 'isKnownVariant', 'ESP_AF', 'TG_AF', 'CCDS', 'GeneName', 'cDNApos', 'relcDNApos', 'CDSpos', 'relCDSpos', 'protPos', 'relProtPos', 'Domain', 'Dst2Splice', 'Dst2SplType', 'Exon', 'Intron', 'oAA', 'nAA', 'Grantham', 'SIFTcat', 'SIFTval', 'ExAC-AF', 'ExAC-ObsAlleles', 'ExAC-TotalAlleles', 'ExAC-HomCount', 'ExAC-HetCount','clinvar','clinvar_simple','rev_comp','pos_nudge','CDSpos_nudge','within_2bp_of_pam_edit'])
        add_to_header_3 = '\tneg_snv_error_ratio\tpre_snv_error_ratio\tpost_snv_error_ratio\texp_snv_error_reads_pre\texp_snv_error_reads_post\tcorr_pre\tcorr_post\tcorr_pre_pseudo\tcorr_post_pseudo\tcorr_pseudo_post_pre_ratio'
        header_string += (add_to_header_1+add_to_header_2+add_to_header_3+'\n')

        #process the variants for SNV edits, and then CADD data if warranted.
        for variant in var_dict:
            skip_clin_var_add = False
            variant_data = var_dict[variant]            
            #get the raw read counts
            read_counts  = {}
            read_counts['pre'] = int(variant_data[index_pre])
            read_counts['post'] = int(variant_data[index_post])
            read_counts['lib'] = int(variant_data[index_lib])
            read_counts['neg'] = int(variant_data[index_neg])
            read_counts['rna'] = int(variant_data[index_rna])
            if variant_data[index_all_hdr_variant] == 'True':
                for timepoint in ['pre','post','lib','neg','rna']:
                    if read_counts[timepoint] > 0:
                        summary_data[exp][timepoint]['oligo_variants'] += 1
            #define requirements for each class of mutation:
            #NOTE - WE're ignoring 'partial' hdr at a single end (because these events are impossible for most exps, anyway)
            c_hdr = (variant_data[index_hdr5] == 'yes' and variant_data[index_hdr3] == 'yes')
            p_hdr_5 = (variant_data[index_hdr5] == 'yes' and variant_data[index_hdr3] == 'no')
            p_hdr_3 = (variant_data[index_hdr5] == 'no' and variant_data[index_hdr3] == 'yes')
            no_hdr = (variant_data[index_hdr5] == 'no' and variant_data[index_hdr3] == 'no')
            no_indels = (int(variant_data[index_del_count]) == 0 and int(variant_data[index_ins_count]) == 0)

            #remember, these snv's are only those NOT in the HDR signature for the experiment AND in the mutated region

            single_snv = (int(variant_data[index_mismatch_count]) == 1)

            if single_snv: #location check on single snv's to prevent editing error's from getting counted here, and trying to get CADD annotations
                snv_edit = get_snv_mismatch(variant_data[index_edit_string],exp_hdr_edits)
                snv_location = int(snv_edit[:snv_edit.find('-X')])
                if exp_mut_start <= snv_location <= exp_mut_end:
                    single_snv = True
                else:
                    single_snv = False

            no_snv = (int(variant_data[index_mismatch_count]) == 0)

            #higher order combinations (mutually exclusive!!):
            wt_wt = (no_hdr and no_snv and no_indels)
            wt_snv = (no_hdr and no_indels and single_snv)
            c_hdr_wt = (c_hdr and no_snv and no_indels)
            c_hdr_snv = (c_hdr and single_snv and no_indels)
            p_hdr_5_wt = (p_hdr_5 and no_snv and no_indels)
            p_hdr_5_snv = (p_hdr_5 and single_snv and no_indels)
            p_hdr_3_wt = (p_hdr_3 and no_snv and no_indels)
            p_hdr_3_snv = (p_hdr_3 and single_snv and no_indels)
            wt_single_del = (int(variant_data[index_del_count]) == 1 and no_hdr)
            c_hdr_single_del = (int(variant_data[index_del_count]) == 1 and c_hdr and (no_snv or single_snv))

            #put into a list
            edit_cat_data = [c_hdr,p_hdr_5,p_hdr_3,no_hdr,no_indels,single_snv,no_snv,wt_wt,wt_snv,c_hdr_wt,c_hdr_snv,p_hdr_5_wt,p_hdr_5_snv,p_hdr_3_wt,p_hdr_3_snv,wt_single_del,c_hdr_single_del]

            var_dict[variant] += edit_cat_data

            # Label SNVs if in mutated region or not (requires more annotation)
            if wt_wt: 
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['wt_wt_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['wt_wt_variants_reads'] += read_counts[timepoint]
                    #Format for extending by CADD data: append a bunch of NA's to the variants entry in the var_dict.
                variant_data += ['REF']*48
                var_dict[variant] = variant_data
                wt_wt_variants.append(variant)

            elif wt_snv: #these guys to calculate your ber base error rates from sequencing (will be slightly inflated by p_hdr_snv's --> wt_snv's)
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['wt_snv_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['wt_snv_reads'] += read_counts[timepoint]
                snv_edit = get_snv_mismatch(variant_data[index_edit_string],exp_hdr_edits) 
                cadd_lookup_key = mismatch_to_cadd_lookup(snv_edit,seq_start,orientation)
                cadd_variant_data = get_cadd_data_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                ref_base = cadd_variant_data[2]
                alt_base = cadd_variant_data[3]
                variant_data += cadd_variant_data
                skip_clin_var_add = False
                if not skip_clin_var_add:
                    if cadd_lookup_key in clinvar_dict:
                        variant_data+=clinvar_dict[cadd_lookup_key]
                    else:
                        variant_data+=['absent','absent']
                variant_data+= get_pos_nudge(ref_base,alt_base)
                within_2_bp = 'False'
                for pam_edit_location in all_pam_edit_locations:
                    if abs(snv_location-pam_edit_location) < 3:
                        within_2_bp = 'True'
                        break
                    else:
                        continue
                variant_data.append(within_2_bp)
                var_dict[variant] = variant_data
                wt_snv_variants.append(variant)

            elif c_hdr_wt: #this gives a count from which to predict the false number of c_hdr_snv's
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['c_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['c_hdr_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['c_hdr_wt_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['c_hdr_wt_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                variant_data += ['REF']*48
                var_dict[variant] = variant_data
                c_hdr_wt_variants.append(variant)

            elif p_hdr_5_wt: #unsure if this needs to be accounted for or not, depends on seq-error
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['p_hdr_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_5_wt_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_5_wt_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                #find the position of the 3' SNV if in the library and change it

                #this line of code will instead add 'NA's
                variant_data += ['REF']*48
                var_dict[variant] = variant_data

            elif p_hdr_3_wt: #probably just for record keeping unless cutting 3,
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['p_hdr_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_3_wt_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_3_wt_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                variant_data += ['REF']*48
                var_dict[variant] = variant_data

            elif c_hdr_snv: #annotate w. cadd
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['c_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['c_hdr_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['c_hdr_snv_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['c_hdr_snv_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                
                snv_edit = get_snv_mismatch(variant_data[index_edit_string],exp_hdr_edits)
                cadd_lookup_key = mismatch_to_cadd_lookup(snv_edit,seq_start,orientation)
                cadd_variant_data = get_cadd_data_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                ####
                #only place this code in c_hdr_snv context for now... should also be applied to correct p_hdr_scenarios if comparing.
                skip_clin_var_add = False
                snv_location = int(snv_edit[:snv_edit.find('-X')])
                CDSpos = get_CDSpos_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                if CDSpos in pam_codon_data:
                    codon_data = pam_codon_data[CDSpos]
                    codon_index = codon_data[0] 
                    wt_codon = codon_data[1]
                    snv_wt_codon = wt_codon[:codon_index]+snv_edit[-1]+wt_codon[codon_index+1:]
                    snv_wt_aa = translate_codon(snv_wt_codon)
                    lib_codon = codon_data[2]
                    snv_lib_codon = lib_codon[:codon_index]+snv_edit[-1]+lib_codon[codon_index+1:]
                    snv_lib_aa = translate_codon(snv_lib_codon)
                    
                    if snv_lib_aa != snv_wt_aa:
                        #can be nonsense or missense then...
                        skip_clin_var_add = True
                        if snv_lib_aa == '*': #what to do if this creates a nonsense codon: 
                            #manuall define the cadd data:
                            non_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],'STOP_GAINED','stop_gained','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],snv_lib_aa,'NA','NA','NA','NA','NA','NA','NA','NA']
                            cadd_variant_data = non_CADD_data+['absent','absent'] #for two clinvar categories

                        else:
                            mis_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],'NON_SYNONYMOUS','missense','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],snv_lib_aa,'NA','NA','NA','NA','NA','NA','NA','NA']
                            cadd_variant_data = mis_CADD_data+['absent','absent']
                            
                else:
                    pass

                ref_base = cadd_variant_data[2]
                alt_base = cadd_variant_data[3]
                variant_data += cadd_variant_data
                if not skip_clin_var_add:
                    if cadd_lookup_key in clinvar_dict:
                        variant_data+=clinvar_dict[cadd_lookup_key]
                    else:
                        variant_data+=['absent','absent']
                variant_data+= get_pos_nudge(ref_base,alt_base)
                within_2_bp = 'False'
                for pam_edit_location in all_pam_edit_locations:
                    if abs(snv_location-pam_edit_location) < 3:
                        within_2_bp = 'True'
                        break
                    else:
                        continue
                variant_data.append(within_2_bp)
                var_dict[variant] = variant_data
                ###

                c_hdr_snv_variants.append(variant)

            elif p_hdr_5_snv: #annotate w. cadd
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['p_hdr_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_5_snv_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_5_snv_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                snv_edit = get_snv_mismatch(variant_data[index_edit_string],exp_hdr_edits) 
                cadd_lookup_key = mismatch_to_cadd_lookup(snv_edit,seq_start,orientation)
                cadd_variant_data = get_cadd_data_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                ####
                #only place this code in c_hdr_snv context for now... should also be applied to correct p_hdr_scenarios if comparing.
                skip_clin_var_add = False
                snv_location = int(snv_edit[:snv_edit.find('-X')])
                CDSpos = get_CDSpos_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                if (CDSpos in pam_codon_data) and (abs(snv_location-edit_5_locations[0])<6): #5' edits must be close to the snv_edit.
                    # check the translation
                    codon_data = pam_codon_data[CDSpos]
                    codon_index = codon_data[0] 
                    wt_codon = codon_data[1]
                    snv_wt_codon = wt_codon[:codon_index]+snv_edit[-1]+wt_codon[codon_index+1:]
                    snv_wt_aa = translate_codon(snv_wt_codon)
                    lib_codon = codon_data[2]
                    snv_lib_codon = lib_codon[:codon_index]+snv_edit[-1]+lib_codon[codon_index+1:]
                    snv_lib_aa = translate_codon(snv_lib_codon)
                    
                    if snv_lib_aa != snv_wt_aa:
                        #can be nonsense or missense then...
                        skip_clin_var_add = True
                        if snv_lib_aa == '*': #what to do if this creates a nonsense codon: 
                            #manuall define the cadd data:
                            non_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],'STOP_GAINED','stop_gained','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],snv_lib_aa,'NA','NA','NA','NA','NA','NA','NA','NA']
                            cadd_variant_data = non_CADD_data+['absent','absent'] #for two clinvar categories
                            
                        else:
                            mis_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],'NON_SYNONYMOUS','missense','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],snv_lib_aa,'NA','NA','NA','NA','NA','NA','NA','NA']
                            cadd_variant_data = mis_CADD_data+['absent','absent']
                            
                else:
                    pass
                ref_base = cadd_variant_data[2]
                alt_base = cadd_variant_data[3]
                variant_data += cadd_variant_data
                if not skip_clin_var_add:
                    if cadd_lookup_key in clinvar_dict:
                        variant_data+=clinvar_dict[cadd_lookup_key]
                    else:
                        variant_data+=['absent','absent']
                variant_data+= get_pos_nudge(ref_base,alt_base)
                within_2_bp = 'False'
                for pam_edit_location in all_pam_edit_locations:
                    if abs(snv_location-pam_edit_location) < 3:
                        within_2_bp = 'True'
                        break
                    else:
                        continue
                variant_data.append(within_2_bp)
                var_dict[variant] = variant_data
                ###


            elif p_hdr_3_snv: #annotate w. cadd
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['p_hdr_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_3_snv_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_3_snv_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]
                snv_edit = get_snv_mismatch(variant_data[index_edit_string],exp_hdr_edits) 
                cadd_lookup_key = mismatch_to_cadd_lookup(snv_edit,seq_start,orientation)
                cadd_variant_data = get_cadd_data_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                ####
                skip_clin_var_add = False
                snv_location = int(snv_edit[:snv_edit.find('-X')])
                CDSpos = get_CDSpos_for_variant(cadd_lookup_key,cadd_data_dod[amplicon])
                if (CDSpos in pam_codon_data) and (abs(snv_location-edit_3_locations[0])<6): #3' edits must be close to the snv_edit.
                    # check the translation
                    codon_data = pam_codon_data[CDSpos]
                    codon_index = codon_data[0] 
                    wt_codon = codon_data[1]
                    snv_wt_codon = wt_codon[:codon_index]+snv_edit[-1]+wt_codon[codon_index+1:]
                    snv_wt_aa = translate_codon(snv_wt_codon)
                    lib_codon = codon_data[2]
                    snv_lib_codon = lib_codon[:codon_index]+snv_edit[-1]+lib_codon[codon_index+1:]
                    snv_lib_aa = translate_codon(snv_lib_codon)
                    
                    if snv_lib_aa != snv_wt_aa:
                        #can be nonsense or missense then...
                        skip_clin_var_add = True
                        if snv_lib_aa == '*': #what to do if this creates a nonsense codon: 
                            #manuall define the cadd data:
                            non_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],'STOP_GAINED','stop_gained','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],snv_lib_aa,'NA','NA','NA','NA','NA','NA','NA','NA']
                            cadd_variant_data = non_CADD_data+['absent','absent'] #for two clinvar categories 
                            print 'ADDED 2 NAs from p_hdr_3_snv'
                            print cadd_variant_data

                        else:
                            mis_CADD_data = [cadd_data[cadd_lookup_key][0],cadd_data[cadd_lookup_key][1],cadd_data[cadd_lookup_key][2],cadd_data[cadd_lookup_key][4],'NON_SYNONYMOUS','missense','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',cadd_data[cadd_lookup_key][94],cadd_data[cadd_lookup_key][95],cadd_data[cadd_lookup_key][96],cadd_data[cadd_lookup_key][97],cadd_data[cadd_lookup_key][98],cadd_data[cadd_lookup_key][99],cadd_data[cadd_lookup_key][100],cadd_data[cadd_lookup_key][101],cadd_data[cadd_lookup_key][102],cadd_data[cadd_lookup_key][103],cadd_data[cadd_lookup_key][104],cadd_data[cadd_lookup_key][105],cadd_data[cadd_lookup_key][106],cadd_data[cadd_lookup_key][107],snv_lib_aa,'NA','NA','NA','NA','NA','NA','NA','NA']
                            cadd_variant_data = mis_CADD_data+['absent','absent']
                            print 'ADDED 2 NAs from p_hdr_3_snv'
                            print cadd_variant_data
                else:
                    pass
                ref_base = cadd_variant_data[2]
                alt_base = cadd_variant_data[3]
                variant_data += cadd_variant_data
                if len(cadd_variant_data) == 39:
                    print cadd_variant_data

                if not skip_clin_var_add:
                    if cadd_lookup_key in clinvar_dict:
                        variant_data+=clinvar_dict[cadd_lookup_key]
                    else:
                        variant_data+=['absent','absent']

                variant_data+= get_pos_nudge(ref_base,alt_base)
                within_2_bp = 'False'
                for pam_edit_location in all_pam_edit_locations:
                    if abs(snv_location-pam_edit_location) < 3:
                        within_2_bp = 'True'
                        break
                    else:
                        continue
                variant_data.append(within_2_bp)
                var_dict[variant] = variant_data
                ###

            elif wt_single_del: #just for nhej reference
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                variant_data += ['NA']*48
                var_dict[variant] = variant_data

            elif c_hdr_single_del: #just for controls
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['c_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['c_hdr_reads'] += read_counts[timepoint]
                variant_data += ['NA']*48
                var_dict[variant] = variant_data

            #BE SURE TO COUNT ALL THE HDR reads even if they aren't perfect SNVs or WT (above)
            elif (p_hdr_5 or p_hdr_3) and not (p_hdr_5_wt or p_hdr_5_snv or p_hdr_3_wt or p_hdr_3_snv):
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['p_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['p_hdr_reads'] += read_counts[timepoint]
                    if no_indels:
                        summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                        summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint]                        
                variant_data += ['NA']*48
                var_dict[variant] = variant_data

            elif c_hdr and not (c_hdr_wt or c_hdr_snv):
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    summary_data[exp][timepoint]['c_hdr_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                    summary_data[exp][timepoint]['c_hdr_reads'] += read_counts[timepoint]
                    if no_indels:
                        summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                        summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint] 
                variant_data += ['NA']*48
                var_dict[variant] = variant_data

            else: #fill in all the stats categories with place holders.
                for timepoint in ['pre','post','lib','neg','rna']:
                    summary_data[exp][timepoint]['total_reads'] += read_counts[timepoint]
                    if no_indels:
                        summary_data[exp][timepoint]['full_length_variants'] += increment_variant_count(var_dict[variant][header_index(header_string,timepoint)])
                        summary_data[exp][timepoint]['full_length_reads'] += read_counts[timepoint] 
                variant_data += ['NA']*48
                var_dict[variant] = variant_data

        ### END OF COUNTING VARIANT TYPES AND ASSIGNING CADD SCORES
        ### CODING OUT LIBRARY CORRECTIONS FOR NOW!
        # Calculate per base error rates, append to all SNV reads
        #expand this to calculate per base error reads for all p_hdr_snv_s as well (unsure how exactly)
        wt_wt_reads_neg = int(summary_data[exp]['neg']['wt_wt_variants_reads'])
        wt_wt_reads_pre = int(summary_data[exp]['pre']['wt_wt_variants_reads'])
        wt_wt_reads_post = int(summary_data[exp]['post']['wt_wt_variants_reads'])
        c_hdr_wt_reads_pre = int(summary_data[exp]['pre']['c_hdr_wt_reads'])
        c_hdr_wt_reads_post = int(summary_data[exp]['post']['c_hdr_wt_reads'])
        c_hdr_wt_reads_lib = int(summary_data[exp]['lib']['c_hdr_wt_reads'])
        c_hdr_wt_reads_rna = int(summary_data[exp]['rna']['c_hdr_wt_reads'])
        for variant in var_dict:
            if variant in c_hdr_snv_variants:
                #get the raw counts for the variant:
                c_hdr_snv_pre = int(var_dict[variant][index_pre])
                c_hdr_snv_post = int(var_dict[variant][index_post])
                c_hdr_snv_rna = int(var_dict[variant][index_rna])
                c_hdr_snv_lib = int(var_dict[variant][index_lib])
                pre_snv_error_ratio = 0.0
                post_snv_error_ratio = 0.0
                lib_snv_error_ratio = 0.0
                rna_snv_error_ratio = 0.0
                #get the snv:
                snv_edit = get_snv_mismatch(var_dict[variant][index_edit_string],exp_hdr_edits)
                snv_edit_split = snv_edit.strip().split('-')
                snv_location = int(snv_edit_split[0])
                snv_base = snv_edit_split[2]
                wt_snv_seq = ref_seq[:snv_location-1]+snv_base+ref_seq[snv_location:]
                wt_snv_pre = 0
                wt_snv_post = 0
                wt_snv_rna = 0
                wt_snv_lib = 0
               # wt_snv_lib = 0
                if wt_snv_seq in wt_snv_variants:
                    wt_snv_pre = int(var_dict[wt_snv_seq][index_pre])
                    wt_snv_post = int(var_dict[wt_snv_seq][index_post])
                    wt_snv_lib = int(var_dict[wt_snv_seq][index_lib])
                    wt_snv_rna = int(var_dict[wt_snv_seq][index_rna])
                    wt_snv_neg = int(var_dict[wt_snv_seq][index_neg])
                else:
                    #print 'SNV edit ' + str(snv_edit) + ' not detected in wt context for error correction in neg. sample.'
                    pass
                if wt_wt_reads_neg == 0:
                    print 'Error! No reference alleles detected in neg sample for error correction.'
                neg_snv_error_ratio = wt_snv_neg/ float(wt_wt_reads_neg)
                pre_snv_error_ratio = wt_snv_pre/ float(wt_wt_reads_pre)
                post_snv_error_ratio = wt_snv_post/ float(wt_wt_reads_post)
                #THESE ARE BASED ON THE NEGATIVE CONTROL SAMPLE NOW!
                exp_snv_error_reads_pre = int(c_hdr_wt_reads_pre * neg_snv_error_ratio)
                exp_snv_error_reads_post = int(c_hdr_wt_reads_post * neg_snv_error_ratio)
                #exp_snv_error_reads_lib = c_hdr_snv_lib * ave_snv_error_ratio

                corr_pre = c_hdr_snv_pre - exp_snv_error_reads_pre
                corr_post = c_hdr_snv_post - exp_snv_error_reads_post
                if corr_pre < 0:
                    corr_pre = 0
                if corr_post < 0:
                    corr_post = 0
                corr_pre_pseudo = corr_pre + 1
                corr_post_pseudo = corr_post + 1
                corr_pseudo_post_pre_ratio = float(corr_post_pseudo)/corr_pre_pseudo
                #
                var_dict[variant].extend([neg_snv_error_ratio,pre_snv_error_ratio,post_snv_error_ratio,exp_snv_error_reads_pre,exp_snv_error_reads_post,corr_pre,corr_post,corr_pre_pseudo,corr_post_pseudo,corr_pseudo_post_pre_ratio])
            else:
                var_dict[variant].extend(['NA']*10)

        output_file_name = i[:-4]+'_annotated.txt'
        output_file = open(my_dir+output_subdir+'/'+output_file_name, 'w')
        output_file.write(header_string)
        for variant in ordered_variants:
            output_data = var_dict[variant]
            output_strings = [variant]
            for k in output_data:
                output_strings.append(str(k))
            output_file.write('\t'.join(output_strings)+'\n')

        summary_header_titles = ['total_reads','c_hdr_variants','c_hdr_reads','c_hdr_wt_variants','c_hdr_wt_reads','c_hdr_snv_variants','c_hdr_snv_reads','p_hdr_variants','p_hdr_reads','t_hdr_variants','t_hdr_reads','p_hdr_5_snv_reads','p_hdr_5_snv_variants','p_hdr_5_wt_reads','p_hdr_5_wt_variants','p_hdr_3_snv_reads','p_hdr_3_snv_variants','p_hdr_3_wt_reads','p_hdr_3_wt_variants','full_length_variants','full_length_reads','wt_wt_variants','wt_wt_variants_reads','wt_snv_variants','wt_snv_reads','wt_wt_freq','full_length_freq','c_hdr_freq','c_hdr_wt_freq','p_hdr_freq','p_hdr_5_wt_freq','p_hdr_3_wt_freq','t_hdr_rate','hq_snv_hdr_rate','oligo_variants','possible_oligo_variants']
        summary_output_file_name = exp+'_summary.txt'
        summary_output_file = open(my_dir+output_subdir+'/'+summary_output_file_name,'w')
        #write the header once
        summary_output_file.write('exp\ttimepoint\t'+'\t'.join(summary_header_titles)+'\n')

        for timepoint in ['pre','post','lib','neg','rna']:
            output_list = [exp,timepoint]
            #Updating the t_hdr stats
            summary_data[exp][timepoint]['t_hdr_variants'] = summary_data[exp][timepoint]['p_hdr_variants'] + summary_data[exp][timepoint]['c_hdr_variants']
            summary_data[exp][timepoint]['t_hdr_reads'] = summary_data[exp][timepoint]['p_hdr_reads'] + summary_data[exp][timepoint]['c_hdr_reads']
            #Calculating the t_hdr_rate and the hq_snv_hdr_rate
            summary_data[exp][timepoint]['t_hdr_rate'] = float(summary_data[exp][timepoint]['t_hdr_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['hq_snv_hdr_rate'] = float(summary_data[exp][timepoint]['c_hdr_snv_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['wt_wt_freq'] = float(summary_data[exp][timepoint]['wt_wt_variants_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['full_length_freq'] = float(summary_data[exp][timepoint]['full_length_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['c_hdr_freq'] = float(summary_data[exp][timepoint]['c_hdr_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['c_hdr_wt_freq'] = float(summary_data[exp][timepoint]['c_hdr_wt_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['p_hdr_freq'] = float(summary_data[exp][timepoint]['p_hdr_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['p_hdr_5_wt_freq'] = float(summary_data[exp][timepoint]['p_hdr_5_wt_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['p_hdr_3_wt_freq'] = float(summary_data[exp][timepoint]['p_hdr_3_wt_reads'])/summary_data[exp][timepoint]['total_reads']
            summary_data[exp][timepoint]['possible_oligo_variants'] = len_oligo_variants
            for header in summary_header_titles:
                output_list.append(str(summary_data[exp][timepoint][header]))
            summary_output_file.write('\t'.join(output_list)+'\n')

        #go line by line in the summary dictionary for the experiment (not the whole sequencing run -- CORRECT) and spit out all the stats in a matrix form
print 'Time elapsed', (datetime.now()-startTime)