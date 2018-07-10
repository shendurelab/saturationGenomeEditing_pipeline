#sam_to_edits_wRNA_pipeline_20170517.py

#NOTE - THIS VERSION UPDATED FROM PREVIOUS to take any variables as sys.argv's
#Updates for RNA version:
	#threshld (1 in 1e6) req. in either RNA, pre, post, or lib.
	#5 sam files per experiment:  pre post lib neg rna
	#adds 21 additional columns to output:  1.) whether or not variant is in the oligo library ('all_hdr_variant','True'/'False')
		#if the variant is a cHDR variant (all PAM edits present), then script searches for pHDR variant (only PAM edit at cut site) by using the reference base to fill in the other HDR site.
		#adds 20 columns for paired pHDR variant (raw reads * 5)(tHDR reads * 5)(tHDR frequencies * 5)(tHDR ratios *5)
		#adds all 'NA's if not cHDR, unless it is in the all_hdr_variants list, then it adds the same measurements for comparison without any adds from pHDR equivalent (because there is none -- i.e. a SNV at PAM edit site...)
import sys
import os
import subprocess
import operator

#amplicon_list = sys.argv[1] ## i.e. X4+X5+X17+... 
#exp_groupings = sys.argv[2] ##

#--------------------------------------------------------------------------------------------------------------------------------------------
#HARCODEDED SECTION
'''	1. Reference file
	2. Editing info file, and info for which exp goes with which editing info.
	3. Read threshold -- currently based on one read per hundred thousand.
	4. Alignment score threshold (below).
	5. Experimental groupings of sam files.
'''
#Harcoded paths and sets of samples:
#SAMPLE SAM files in:
my_dir =  os.getcwd()
#could switch to populate this list from the file names in the sam folder, too. Hardcode for now.
amplicon_list = sys.argv[1].split('+')

#REFERENCE FILES (in fasta format -- takes second line of file), hardcoded to look in main BRCA1 fasta storage:
fasta_dir = '/net/shendure/vol10/projects/SGE/nobackup/BRCA1/fasta'
ref_seqs = {}
for amplicon in amplicon_list:
	ref_file = open(fasta_dir+'/'+amplicon+'.fa', 'r')
	ref_header = ref_file.readline()
	ref_seq = ref_file.readline().strip()
	ref_seqs[amplicon] = ref_seq
	ref_file.close()

#This can stay hardcoded for all BRCA1 experiments, now.
edits_per_amp = {}
editing_info = open('/net/shendure/vol10/projects/SGE/nobackup/BRCA1/BRCA1_editing_data.txt', 'r')
editing_info_header = editing_info.readline()
for line in editing_info:
	edits = line.strip().split()
	#each amplicon in the file is a key that points to a list containing the info for that amplicon [5' HDR, 3' HDR, cut_pos, cut_end, mut_start, mut_end] 
	edits_per_amp[edits[0]] = edits[1:]


#sets a reads threshold for making it into the output
read_threshold = .000002 #reads to report a variant observed -- set to a fraction (such as .000001 for 1 in 1M)
#ALIGNMENT SCORE THRESHOLD:
alignment_score_threshold = 300
#define experimental sets: format is baseline sample-->[list of corresponding samples ordered pre, post, lib, neg ]

exp_groupings_list = sys.argv[2].split('+')
exp_groupings = {}
for grouping in exp_groupings_list:
	g_info = grouping.split(',')
	exp_groupings[g_info[0]] = g_info[1:] 

print exp_groupings
#DEFINE FUNCTIONS:

#--------------------------------------------------------------------------------------------------------------------------------------------
#cigar_to_edits -- A function that takes in a cigar string, a variant sequence, and a reference, and returns 'expanded' cigar and a list of edits
#call it on each variant after pile-up
def cigar_to_edits(cigar, seq, ref_seq):
	MDI_indexes = []
	for x in range(0,len(cigar)):
		if (cigar[x] == "M") or (cigar[x] == "D") or (cigar[x] == "I"):
			MDI_indexes.append(x)
	#now iterate over the MDI index list and do something different to call variants for each of the different options
	#these idexes will keep track of where we are in the read ('seq'), in the ref_seq, and the cigar
	cigar_index = 0
	read_index = 0
	ref_index = 0
	#keeps matches listed
	mod_cigar = ''
	#only shows location of edits
	edit_string = ''
	for mdi_index in MDI_indexes:
		char = cigar[mdi_index]
		length = int(cigar[cigar_index:mdi_index])
		#cigar_index no longer useable for this chunk after this
		cigar_index = mdi_index+1
		#will call each SNV individually
		#X will be used to indicate a mismatch (an isolated X is a SNV - not specified here, though)
		# mod_cigar FORMAT: [bases][type][if insertion or mismatch, base id(s)]
		# edit_string FORMAT: comma-separated list of variants from WT:  [position]-[variant type]-[length (for I/D)]-[base id(s) (for I/X)]
		if char == 'M':
			match_length = 0
			#x will be relative to the reference
			#need an adjustment factor to get the appropriate index in the seq
			adj = ref_index-read_index
			for x in range(ref_index,ref_index+length):
				if seq[x-adj] == ref_seq[x]:
					match_length+=1
					if x == (ref_index+length-1):
						mod_cigar += (str(match_length)+'M')
						match_length = 0
				elif seq[x-adj] != ref_seq[x]:
					if match_length >= 1:
						mod_cigar += (str(match_length)+'M'+'1X'+seq[x-adj])
					elif match_length == 0:
						mod_cigar += ('1X'+seq[x-adj])
					edit_string += (str(ref_index+x+1)+'-'+'X-'+seq[x-adj]+',')
					match_length = 0
			#adjust the read_index and ref_index
			ref_index += length
			read_index += length
			
		elif char == 'D':
			mod_cigar += (str(length)+'D')
			edit_string += (str(ref_index+1)+'-D'+str(length)+',')
			#call position of the deleted bases
			#adjust the ref_index
			ref_index += length
		elif char == 'I':
			insertion_sequence = seq[read_index:read_index+length]
			mod_cigar += (str(length)+'I'+insertion_sequence)
			edit_string += (str(ref_index+1)+'-I'+str(length)+'-'+insertion_sequence+',')
			#adjust the read_index (ref_index stays the same)
			read_index += length
	if edit_string == '':
		edit_string = '-WT'
	else: #strip off the last comma
		edit_string = edit_string[:-1]
	#store the mod_cigar and the edit_string in the variant_dict
	#also return the number of edits, and the number of each type of edit
	return([mod_cigar, edit_string])

#--------------------------------------------------------------------------------------------------------------------------------------------

#dol_lookup -- a function to lookup a key-number pairing in a dicitonary -- instead of returning a key error, return a 0 if not present)
def dol_lookup(counts_dol,key): #will return 0 for the value if not present, otherwise the first item for each list with v:[counts,others] format
	if key in counts_dol:
		return counts_dol[key][0]
	else:
		return 0

#--------------------------------------------------------------------------------------------------------------------------------------------

#returns a dictionary for all variants, formatted {variant:[counts,cigar]} -- demands all variants have alignment score over min
def sam_to_variant_cigar_dict(sam_file, min_alignment_score):
	line_count = -2
	variant_dict = {}
	reads_not_aligning = 0
	for line in sam_file:
		if line_count < 0:
			line_count += 1
			continue
		else:
			line_count += 1
			sam_data = line.strip().split('\t')
			#cigar string
			cigar = sam_data[5]
			#this is the alignment score -- not necessarily useful unless mispriming common
			AS = sam_data[11]
			#sequence of read
			seq = sam_data[9]
			#filter based on the alignment score -- good alignments should be over ~500
			index_of_score = AS.find('i:')+2
			#this is the actual score
			score = float(AS[index_of_score:])
			#hard coded alignment score cuttoff (using needleall with 10 gap open and 0.5 gap extension -- 100 is fairly arbitrary -- very low to be inclusive. could raise to only look at better reads)
			if score > min_alignment_score:
				if seq in variant_dict:
					variant_dict[seq][0]+=1
				else:
					variant_dict[seq]=[1,cigar]
			else: 
				reads_not_aligning += 1
	print reads_not_aligning, 'reads not aligning in sam file "', str(sam_file), ' using min_align_score of ', min_alignment_score, ' out of ', line_count, ' total reads.' 
	return variant_dict #returns the final dictionary after parsing
#-------------------------------------------------------------------------------------------------------------------------------------------


def sort_counts_dol(counts_dol): #specifies count dict is variant --> [counts, as part of a larger list]
    return sorted(counts_dol.keys(), key=lambda k: -1*int(counts_dol[k][0]))

#--------------------------------------------------------------------------------------------------------------------------------------------

def get_read_count_dol(counts_dol): #specifies count dict is variant --> [counts, as part of a larger list]
    total_reads = 0
    for my_key in counts_dol:
        total_reads += counts_dol[my_key][0]
    return total_reads

#--------------------------------------------------------------------------------------------------------------------------------------------

def header_index(header_line, category):
	header_list = header_line.strip().split('\t') #header list will have 'variant' category as well
	var_dict_index = header_list.index(category)-1
	return var_dict_index


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

#--------------------------------------------------------------------------------------------------------------------------------------------

#to be run in directory with sam files --  must have reference files (.fa) present there, too
working_dir = os.getcwd()
dict_of_var_dicts = {}
sample_list = []

#run on all sam files in the cwd
for i in os.listdir(os.getcwd()):
	if not i.endswith(".sam"):
		continue
	elif i.endswith(".sam"):
		#gets the sample name (i.e X2_HDRL), the timepoint (i.e. pre), and the editing info expected
		index_of_first_dot = i.find('.')
		sample = i[:index_of_first_dot]
		index_of_first_under = i.find('_')
		timepoint = i[index_of_first_under+1:index_of_first_dot]
		experiment = i[:index_of_first_under]
		if amplicon not in amplicon_list:
			print '! Amplicon "'+amplicon+' not detected in amplicon list.'
		if sample not in sample_list:
			sample_list.append(sample)
		#opens the sam_file and imports the count dict:
		with open(i, 'r') as sam_file:
			dict_of_var_dicts[sample] = sam_to_variant_cigar_dict(sam_file,alignment_score_threshold)

'''

'''
#remember, want things that score highly in D11 as well -- would filter on presence either in pre or post
#running header for each file:
header_string = 'variant\tcigar\text_cigar\tedit_string\tpre\tpost\tlib\tneg\trna\tpre_pseudo\tpost_pseudo\tlib_pseudo\tneg_pseudo\trna_pseudo\tpre_freq\tpost_freq\tlib_freq\tneg_freq\trna_freq\tpre_pseudo_freq\tpost_pseudo_freq\tlib_pseudo_freq\tneg_pseudo_freq\trna_pseudo_freq\tpre_lib_ratio\tpost_pre_ratio\tpost_lib_ratio\trna_pre_ratio\trna_post_ratio\t'
header_string+='hdr5\thdr3\t'
header_string+='mismatch_count\tdel_count\tins_count\t'
header_string+='all_hdr_variant\t'
header_string+='pHDR_pre\tpHDR_post\tpHDR_lib\tpHDR_neg\tpHDR_rna\ttHDR_pre\ttHDR_post\ttHDR_lib\ttHDR_neg\ttHDR_rna\t'
header_string+='tHDR_pre_pseudo_freq\ttHDR_post_pseudo_freq\ttHDR_lib_pseudo_freq\ttHDR_neg_pseudo_freq\ttHDR_rna_pseudo_freq\ttHDR_pre_lib_ratio\ttHDR_post_pre_ratio\ttHDR_post_lib_ratio\ttHDR_rna_pre_ratio\ttHDR_rna_post_ratio\t'

for exp in exp_groupings.keys():
	if 'r' in exp:
		amplicon = exp[:exp.find('r')]
	else:
		amplicon = exp
	ref_seq = ref_seqs[amplicon]
	#enumerate the list of all ordered variants:
	#edits per amplicon dict format --> amp:[5' HDR, 3' HDR, cut_pos, cut_end, mut_start, mut_end]

	cSNV_pSNV_pairs = {} # a dict in format cSNV --> pSNV
	#sets the editing info to evaluate the dictionary
	#my_editing_info = edits_per_amp[exp_amplicon_pairings[exp]]
	my_editing_info = edits_per_amp[amplicon]
	#lists of 5' and 3' edits
	edits5 = my_editing_info[0].strip().split(',')
	edits3 = my_editing_info[1].strip().split(',')
	cut_pos = int(my_editing_info[2])
	cut_end = my_editing_info[3]
	mut_start = my_editing_info[4]
	mut_end = my_editing_info[5]
	all_hdr_variants = get_all_hdr_variants(amplicon,ref_seqs,edits_per_amp)
	#relevant count dictionaries for each
	pre_var_dict = dict_of_var_dicts[exp_groupings[exp][0]]
	post_var_dict = dict_of_var_dicts[exp_groupings[exp][1]]
	lib_var_dict = dict_of_var_dicts[exp_groupings[exp][2]]
	neg_var_dict = dict_of_var_dicts[exp_groupings[exp][3]]
	rna_var_dict = dict_of_var_dicts[exp_groupings[exp][4]]

	pre_var_rc = float(get_read_count_dol(pre_var_dict))
	post_var_rc = float(get_read_count_dol(post_var_dict))
	lib_var_rc = float(get_read_count_dol(lib_var_dict))
	neg_var_rc = float(get_read_count_dol(neg_var_dict))
	rna_var_rc = float(get_read_count_dol(rna_var_dict))

	#create a list of all variants present in either pre_var_dict or post_var_dict (#edited here to include lib as well) -- sorted by counts
	pre_variants = sort_counts_dol(pre_var_dict)
	print len(pre_variants), 'pre-variants'
	post_variants = sort_counts_dol(post_var_dict)
	print len(post_variants), 'post-variants'
	lib_variants = sort_counts_dol(lib_var_dict)
	print len(lib_variants), 'lib-variants'
	rna_variants = sort_counts_dol(rna_var_dict)
	print len(rna_variants), 'rna-variants'
	#print post_variants
	#what's going on here?! -- Answer from Aaron;  need to use a dictionary not a list for comparison -- will look up way faster due to indexed nature.
	unique_post_variants = []
	unique_rna_variants = []
	for post_variant in post_variants:
		if post_variant not in pre_var_dict:
			unique_post_variants.append(post_variant)
	for rna_variant in rna_variants:
		if rna_variant not in pre_var_dict:
			unique_rna_variants.append(rna_variant)

	print len(unique_post_variants), 'unique post variants (not in pre)'
	print len(unique_rna_variants), 'unique rna variants (not in pre)'
	pre_post_variants = pre_variants+unique_post_variants
	unique_lib_variants = []
	for lib_variant in lib_variants:
		if (lib_variant not in pre_var_dict) and (lib_variant not in post_var_dict):
			unique_lib_variants.append(lib_variant)
	print len(unique_lib_variants), 'unique lib variants'
	exp_variants = pre_post_variants + unique_lib_variants
	print len(exp_variants), 'total exp. variants (iclu. lib)'

	pre_var_pseudo_rc = pre_var_rc + len(pre_var_dict)
	post_var_pseudo_rc = post_var_rc + len(post_var_dict)
	lib_var_pseudo_rc = lib_var_rc + len(lib_var_dict)
	neg_var_pseudo_rc = neg_var_rc + len(neg_var_dict)
	rna_var_pseudo_rc = rna_var_rc + len(rna_var_dict)

	#place to store aggregate data from each sample:
	exp_master_dict = {}
	pHDR_print_count = 0
	for variant in exp_variants:
		outline_data = []
		if variant in pre_var_dict:
			cigar = pre_var_dict[variant][1]
		elif variant in post_var_dict:
			cigar = post_var_dict[variant][1]
		else:
			cigar = lib_var_dict[variant][1]
		outline_data+=[cigar]
		outline_data+=cigar_to_edits(cigar,variant,ref_seq) #define ref seq here first!!!
		raw_reads = [dol_lookup(pre_var_dict,variant),dol_lookup(post_var_dict,variant),dol_lookup(lib_var_dict,variant),dol_lookup(neg_var_dict,variant),dol_lookup(rna_var_dict,variant)]
		pseudo_reads = []
		for i in raw_reads:
			pseudo_reads.append(1+i)
		raw_freq = [raw_reads[0]/pre_var_rc,raw_reads[1]/post_var_rc,raw_reads[2]/lib_var_rc,raw_reads[3]/neg_var_rc,raw_reads[4]/rna_var_rc]
		pseudo_freq = [pseudo_reads[0]/pre_var_pseudo_rc,pseudo_reads[1]/post_var_pseudo_rc,pseudo_reads[2]/lib_var_pseudo_rc,pseudo_reads[3]/neg_var_pseudo_rc,pseudo_reads[4]/rna_var_pseudo_rc]
		#only 3 ratios calculated, all from pseudocounts
		ratios = [pseudo_freq[0]/pseudo_freq[2],pseudo_freq[1]/pseudo_freq[0],pseudo_freq[1]/pseudo_freq[2],pseudo_freq[4]/pseudo_freq[0],pseudo_freq[4]/pseudo_freq[1]]
		outline_data = outline_data + raw_reads + pseudo_reads + raw_freq + pseudo_freq + ratios

		#test to see if the variants' edits match the HDR edits (from outside tdl file specified above)
		var_edits = outline_data[2].split(',')
		edited_5_count = 0
		edited_3_count = 0

		for edit in edits5:
			if edit in var_edits:
				edited_5_count +=1
		if edited_5_count == len(edits5):
			edited_5 = 'yes'
		elif 0 < edited_5_count < len(edits5):
			edited_5 = 'partial'
		elif edited_5_count == 0:
			edited_5 = 'no'
		for edit in edits3:
			if edit in var_edits:
				edited_3_count+=1
		if edited_3_count == len(edits3):
			edited_3 = 'yes'
		elif 0 < edited_3_count < len(edits3):
			edited_3 = 'partial' 
		elif edited_3_count == 0:
			edited_3 = 'no'
		hdr_data = []
		if cut_end == '5':
			hdr_data = [edited_5,edited_3] #determines the order of the hdr_data list to be [proximal end, distal end] to cut
		elif cut_end == '3':
			hdr_data = [edited_3,edited_5]
		else:
			print cut_end, 'cut_end'
			print 'error:  hdr_data not defined'
		if hdr_data == ['yes','yes']:
			cHDR = True
			pHDR = False
		elif hdr_data == ['yes','no']:
			pHDR = True
			cHDR = False
		else:
			pHDR = False
			cHDR = False
		outline_data += hdr_data

		#ignoring all edits that are HDR markers
		mismatch_count = 0
		del_count = 0
		ins_count = 0

		#check remaining edits in the string to provide additional annotations on the variant's status	X
		#are any of them indels? how many of each?	X
		#is there only 1 substitution, or multiple or none  X
		#is that substitution within the edited margins / exon margins --- ADD THIS LATER

		for edit in var_edits:
			if edit in edits5:
				continue
			elif edit in edits3:
				continue
			else:
				edit_type = edit[edit.index('-')+1]
				if edit_type == 'D':
					del_count += 1
				elif edit_type == 'I':
					ins_count +=1
				elif edit_type == 'X':
					mismatch_count +=1
				elif edit_type == 'W':
					pass	

		outline_data += [mismatch_count, del_count, ins_count]
		if variant in all_hdr_variants:
			outline_data.append('True')
		else:
			outline_data.append('False')			
		if cHDR and (del_count == 0) and (ins_count == 0) and (mismatch_count<len(edits5+edits3)+2): #only cHDR w/ 1 max mismatch
			#find the pSNV version of the variant and get the counts for that w/ dol_lookup
			if cut_end == '5':
				distal_edits = edits3
			elif cut_end == '3':
				distal_edits = edits5
			pHDR_eq_var = variant
			for edit in distal_edits:
				edit_location = int(edit[:edit.find('-')])
				ref_base = ref_seq[edit_location-1]
				pHDR_eq_var = variant[:edit_location-1]+ref_base+variant[edit_location:]
			if pHDR_eq_var == variant:
				print variant, '\n', pHDR_eq_var
				print cigar_to_edits(cigar,variant,ref_seq)
			if pHDR_print_count <= 5:
				print variant
				print pHDR_eq_var
				pHDR_print_count +=1
			pHDR_raw_reads = [dol_lookup(pre_var_dict,pHDR_eq_var),dol_lookup(post_var_dict,pHDR_eq_var),dol_lookup(lib_var_dict,pHDR_eq_var),dol_lookup(neg_var_dict,pHDR_eq_var),dol_lookup(rna_var_dict,pHDR_eq_var)]
			tHDR_pseudo_reads = [pseudo_reads[0]+pHDR_raw_reads[0],pseudo_reads[1]+pHDR_raw_reads[1],pseudo_reads[2]+pHDR_raw_reads[2],pseudo_reads[3]+pHDR_raw_reads[3],pseudo_reads[4]+pHDR_raw_reads[4]]
			tHDR_pseudo_freqs = [tHDR_pseudo_reads[0]/pre_var_pseudo_rc,tHDR_pseudo_reads[1]/post_var_pseudo_rc,tHDR_pseudo_reads[2]/lib_var_pseudo_rc,tHDR_pseudo_reads[3]/neg_var_pseudo_rc,tHDR_pseudo_reads[4]/rna_var_pseudo_rc]
			tHDR_ratios = [tHDR_pseudo_freqs[0]/tHDR_pseudo_freqs[2],tHDR_pseudo_freqs[1]/tHDR_pseudo_freqs[0],tHDR_pseudo_freqs[1]/tHDR_pseudo_freqs[2],tHDR_pseudo_freqs[4]/tHDR_pseudo_freqs[0],tHDR_pseudo_freqs[4]/tHDR_pseudo_freqs[1]]
			outline_data = outline_data + pHDR_raw_reads + tHDR_pseudo_reads + tHDR_pseudo_freqs + tHDR_ratios
		elif (variant in all_hdr_variants) and not cHDR: #this will carry over the numbers from above for the all_hdr_variants that aren't cHDR
		#error -- this is carrying numbers over for pHDR variants that will not show up at all in RNA because they don't have a PAM edit in intron.
			#determine if remaining pSNV edit is in the mutated region or adapter region (if adapter -- change it back to WT to get RNA reads out)
			if cut_end == '5':
				distal_edits = edits3
			elif cut_end == '3':
				distal_edits = edits5
			pHDR_RNA_var = variant
			for edit in edits3+edits5:
				edit_location = int(edit[:edit.find('-')])
				if not (int(mut_start) <= edit_location <= int(mut_end)):
					ref_base = ref_seq[edit_location-1]
					pHDR_RNA_var = pHDR_RNA_var[:edit_location-1]+ref_base+pHDR_RNA_var[edit_location:] #substitutes in the ref. base at any HDR edit outside of SGE region
			if pHDR_RNA_var != variant:
				print variant, 'variant\n', pHDR_RNA_var, 'pHDR_RNA_var'
				pHDR_raw_reads = [dol_lookup(pre_var_dict,pHDR_RNA_var),dol_lookup(post_var_dict,pHDR_RNA_var),dol_lookup(lib_var_dict,pHDR_RNA_var),dol_lookup(neg_var_dict,pHDR_RNA_var),dol_lookup(rna_var_dict,pHDR_RNA_var)]
				tHDR_pseudo_reads = [pseudo_reads[0]+pHDR_raw_reads[0],pseudo_reads[1]+pHDR_raw_reads[1],pseudo_reads[2]+pHDR_raw_reads[2],pseudo_reads[3]+pHDR_raw_reads[3],pseudo_reads[4]+pHDR_raw_reads[4]]
				tHDR_pseudo_freqs = [tHDR_pseudo_reads[0]/pre_var_pseudo_rc,tHDR_pseudo_reads[1]/post_var_pseudo_rc,tHDR_pseudo_reads[2]/lib_var_pseudo_rc,tHDR_pseudo_reads[3]/neg_var_pseudo_rc,tHDR_pseudo_reads[4]/rna_var_pseudo_rc]
				tHDR_ratios = [tHDR_pseudo_freqs[0]/tHDR_pseudo_freqs[2],tHDR_pseudo_freqs[1]/tHDR_pseudo_freqs[0],tHDR_pseudo_freqs[1]/tHDR_pseudo_freqs[2],tHDR_pseudo_freqs[4]/tHDR_pseudo_freqs[0],tHDR_pseudo_freqs[4]/tHDR_pseudo_freqs[1]]
			else:
				pHDR_raw_reads = raw_reads
				tHDR_pseudo_reads = pseudo_reads
				tHDR_pseudo_freqs = pseudo_freq
				tHDR_ratios = ratios
			outline_data = outline_data + pHDR_raw_reads + tHDR_pseudo_reads + tHDR_pseudo_freqs + tHDR_ratios
		else:
			outline_data += ['NA']*5+['NA']*5+['NA']*5+['NA']*5
		exp_master_dict[variant] = outline_data

	outfile_path = my_dir+'/variant_counts_no_thresh/'+exp+'.txt'
	outfile_open = open(outfile_path,'w')
	outfile_open.write(header_string+'\n')
	pre_threshold = read_threshold*pre_var_rc
	post_threshold = read_threshold*post_var_rc
	lib_threshold = read_threshold*lib_var_rc
	rna_threshold = read_threshold*rna_var_rc
	variants_above_threshold = 0
	for variant in exp_variants:
		if (dol_lookup(pre_var_dict,variant)>pre_threshold) or (dol_lookup(post_var_dict,variant)>post_threshold) or (dol_lookup(lib_var_dict,variant)>lib_threshold) or (dol_lookup(rna_var_dict,variant)>rna_threshold):
			variants_above_threshold +=1
			output_list = [variant]+[str(i) for i in exp_master_dict[variant]]
			outfile_open.write('\t'.join(output_list)+'\n')
		else:
			pass
	print 'Finished creating output file for '+exp+' experiment.'
	print variants_above_threshold, "variants above threshold in pre, post, lib, or RNA sample."


