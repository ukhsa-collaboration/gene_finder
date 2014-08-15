#!/usr/bin/env python
import logging
from subprocess import call
import os,sys,inspect,re,subprocess
import logging
from subprocess import call
import os,sys,inspect,re,subprocess
from itertools import groupby, count
from operator import itemgetter
from collections import Counter
from StringIO import StringIO
from lxml import etree
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import pysam


module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
	module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
	if module_folder not in sys.path:
		sys.path.insert(1, module_folder)
import utility_functions
import log_writer



"""
Function
read_mipelup is the main function, calles two function
a. sequence_lengths_for_ref_alleles: extarct allele size
b. parsing_mpileup: extract static coverage for each allelle

The option of the method
fasta_file[str]:
cut_off[str]:
minimum_coverage[str]:
outdir[str]:
workflow_name[str]:
version[str]:
prefix[str]:
logger[str]:

returns(hash)
"""
def read_mipelup(fasta_file,cut_off,minimum_coverage,outdir,workflow_name,version,prefix,logger):

    fai_file = fasta_file + ".fai"
    pileup_file = outdir + "/tmp/" + prefix + ".pileup"
    sorted_bam = outdir + "/tmp/" + prefix + ".sorted.bam"
    size = sequence_lengths_for_ref_alleles(fai_file) #size {'dfr_G': 498}
    mpileup_dictionary = parsing_mpileup(fasta_file,pileup_file,cut_off,minimum_coverage,size,sorted_bam,logger)#mpileup_dictionary {'dfr_G': {'positions_infos': {1: 55, 2: 55,,,,,},
    #positions_low_coverage': {}, 'position_deletions': {}, 'allele_length': 498, 'inserted_nuc': [], '#position_insertions': {}, 'positions_accepted_depth': {1: 55, 2: 55, 3: 56, 4: 55, 5:
    #'position_mismatchs': {}, 'positions_indels_probabilities': {}, 'sequence_raw': {1: 'A', 2: 'T', 3: 'G', 4: 'A', 5:
    return mpileup_dictionary


"""
1. Extract allele size
The option of the method
fai_file[str]:

return

"""
def sequence_lengths_for_ref_alleles(fai_file):
    ## Extract length of ref sequences by parsing SAMtools.fai file and return dictionary={key = string 'gene_name'; value = integer 'sequence length'}
    size = {}
    with open(fai_file) as fai:
        for line in fai:
            fields = line.split('\t')
            size[fields[0]] = int(fields[1])
    return size

def is_integer(element):
    ## check if object is integer
    try:
        i = int(element)
        return True
    except ValueError:
        return False



"""
2.Calls the following functions

a. pileup_extract_information
b.update_infos_with_deletions
c.update_infos_with_coverage
d.complete_missing_seq
e. updated_infos_with_insertions

The option of the method
fasta_file[str]:
pileup_file[str]:
cut_off[str]:
minimum_coverage[str]:
size[str]:
sorted_bam[str]:
logger[str]:

return

"""
def parsing_mpileup(fasta_file,pileup_file,cut_off,minimum_coverage,size,sorted_bam,logger):
    mpileup_info_hash = {}
    with open(pileup_file) as pileup:
        ## split lines in pileup_file by whitespace: group mpileup info by first column info 'allele_name' #=#
        pileup_split = ( x.split() for x in pileup )
        for allele, lines in groupby(pileup_split, itemgetter(0)):#allele dfr_G,  lines # <itertools._grouper object at 0x1767950>
            allele_length = size[allele]#allele_length 498
            ## start description of dictionaries used to collect infos ###
            position_deletions = {}             # keys(positions) ; values(string = deletion patterns) # 
            position_insertions = {}            # keys(positions) ; values(string = insertion patterns) # 
            position_mismatchs = {}             # keys(positions) ; values(list = (ref_nuc and modified_nuc) #
            positions_infos = {}                # keys(positions) ; values(integer = number of reads at positions) #
            positions_indels_probabilities = {} # keys(positions) ; values(float = ratio start_end of read) #
            sequence_raw = {}                   # keys(positions) ; values(string = predicted_nucleotides) #
            positions_mix = {}                  # keys(positions) ; values(dictionary 2 items key=nuc_type value_1=corresponding_ratio)
            positions_low_coverage = {}         # keys(positions) ; values('$')
            ## end of description ###
            ## extract infos by nucleotide position and update dictionaries with infos
            for line in lines:
                details = list(line) #details ['dfr_G', '2', 'T', '55', ',.,.,..,,.,.,,.,.,...,.,..,..,,,,,.,...,..,.,.,.,.,,.,.', 'IFFIFFIFFFBFIFFFIFFFIBFFFIFFIFFFBIIFFIIFFIFFIIFIIFFIFFF']
                ## ensure that each line contain 6 info columns (allele;position;nuc_ref;total_reads;reads_info;reads_quality): this may change with version of samtools
                if len(details) == 6:
                    allele_name,position,nuc_ref,total_reads,reads_info,reads_quality = details
		    #allele_name:dfr_G, #position: 2, #nuc_ref: T, total_reads: 55, #reads_info:,.,.,..,,.,.,,.,.,...,.,..,..,,,,,.,, reads_quality;IFFIFFIFFFBFI
                    ## pileup_extract_information method extract infos from mpileup line by line
                    mapping_type, nuc, nb_reads,ins_sep,insert_pattern,insert_nb,del_sep,deletion_pattern,del_nb,big_indels_indicator,position_mix_det = pileup_extract_information(nuc_ref.upper(), reads_info, reads_quality, allele, allele_length, int(position),total_reads,cut_off,minimum_coverage)
                    #print mapping_type, nuc, nb_reads,ins_sep,insert_pattern,insert_nb,del_sep,deletion_pattern,del_nb,big_indels_indicator,position_mix_det
		    #mapping_type:'match', #nuc:'A', #nb_reads; 55, #ins_sep,:ins, #insert_pattern: 0, #insert_nb:0, #del_sep: del
		    #deletion_pattern:0, #del_nb:0, #big_indels_indicator:100.0, #position_mix_det: {},  #'match', 'A', 55, 'ins', 0, 0, 'del', 0, 0, 100.0, {}
		   
                    if len(position_mix_det) > 0:
                        positions_mix[int(position)] = position_mix_det                
                    if mapping_type == 'match' or mapping_type == 'mismatch' or mapping_type == 'mix' or mapping_type == 'deletion':
                        positions_infos[int(position)] = int(nb_reads)
                    if  mapping_type == 'low_coverage':
                        positions_low_coverage[int(position)] = '$'
                    if mapping_type == 'mismatch' or mapping_type == 'mix': #????? this would have been already picked up earlier, the second if statment 
                        position_mismatchs[int(position)] = [nuc_ref,nuc]
                    if insert_pattern != 0:                        
                        position_insertions[int(position)] = insert_pattern
                    if deletion_pattern != 0:
                        position_deletions[int(position)] = deletion_pattern
                    if int(big_indels_indicator) >= 70:
                        if int(position) != 1 and int(position) != allele_length:
                            positions_indels_probabilities[int(position)] = str(nb_reads) + ":" + str(big_indels_indicator)        
                    sequence_raw[int(position)] = nuc
                else:
                    log_writer.info_header(logger,"SAMtools version 0.1.18 need to be used")#print "SAMtools version 0.1.18 need to be used"
                    sys.exit(1)
            ### refine mpileup infos ###
            if len(positions_infos) != 0:
                max_depth =  max(positions_infos.values())
                ## define threshold for minimum nb of reads acceptable for coverage depth using minimum_coverage value: default=5
                ## minimum_coverage value = percentage and can be changed using argument -m (-min_cov)
                accepted_reads_depth = int(round(float(max_depth*int(minimum_coverage)/100),0))
                ## remove all positions_infos if number of reads is less than accepted_reads_depth:
                positions_with_accepted_depth = dict([(x,positions_infos[x]) for x in positions_infos if int(positions_infos[x]) > accepted_reads_depth])
                total_positions = range(1,allele_length +1)
                info_per_allele = {}      
                if len(positions_with_accepted_depth) != 0:
                    ## update_with deletion info
                    position_infos_up,seq_raw_u1 = update_infos_with_deletions(positions_infos,sequence_raw,position_deletions,positions_with_accepted_depth,total_positions) #???               
                    ## flag positions with no coverage in sequence_raw
                    seq_raw_u2 = update_infos_with_coverage(seq_raw_u1,total_positions,positions_with_accepted_depth)                    
                    
		    ########################
		    # rediat: uses pysam to add sequence at the begining and end
		    ## add sequences at begining/end if info for these regions are missing in mpileup
                    positions_with_accepted_depth_up,position_mismatchs_up,seq_raw_3,positions_indels_probabilities_up = complete_missing_seq(fasta_file,positions_with_accepted_depth,position_mismatchs,positions_indels_probabilities,allele,allele_length,seq_raw_u2,sorted_bam)                   
		    ## update_with deletion info
                    sequence_raw_up,total_inserted_nuc = updated_infos_with_insertions(positions_with_accepted_depth_up,position_insertions,seq_raw_3)
                    ## flag positions with no coverage in positions_with_accepted_depth_up
                    for tot_pos in total_positions: 
                        if positions_with_accepted_depth_up.has_key(tot_pos):
                            pass
                        else:
                            positions_with_accepted_depth_up[int(tot_pos)] = '$'
                            
                    ## start description of dictionaries to be reported for each allele ###    
                    info_per_allele['sequence_raw'] = sequence_raw_up
                    info_per_allele['positions_infos'] = position_infos_up
                    info_per_allele['positions_indels_probabilities'] = positions_indels_probabilities_up
                    info_per_allele['position_deletions'] = position_deletions
                    info_per_allele['position_insertions'] = position_insertions
                    info_per_allele['position_mismatchs'] = position_mismatchs_up            
                    info_per_allele['positions_mix'] = positions_mix
                    info_per_allele['positions_low_coverage'] = positions_low_coverage
                    info_per_allele['allele_length'] = allele_length
                    info_per_allele['positions_accepted_depth'] = positions_with_accepted_depth_up
                    info_per_allele['inserted_nuc'] = total_inserted_nuc
		    
		    #info_per_allele {'positions_infos': {1: 55, 2: 55, 3: 56, 4: 55 }, 'positions_low_coverage': {}, 'position_deletions': {}, 'allele_length': 498, 'inserted_nuc': []
                    ## end of description
                mpileup_info_hash[allele] = info_per_allele # generate a hash and the key is the allele and value is a coverage stat mertics for each allele
		#mpileup_dictionary {'dfr_G': {'positions_infos': {1: 55, 2: 55,,,,,},
		#positions_low_coverage': {}, 'position_deletions': {}, 'allele_length': 498, 'inserted_nuc': [], '
		#position_insertions': {}, 'positions_accepted_depth': {1: 55, 2: 55, 3: 56, 4: 55, 5:
		#'position_mismatchs': {}, 'positions_indels_probabilities': {}, 'sequence_raw': {1: 'A', 2: 'T', 3: 'G', 4: 'A', 5:
		
    return mpileup_info_hash



"""
2.1
Called within the parsing_mpileup function

calles:
1. combine_consecutive_integers_in_list
2. combine_indels_pattern
3. extract_most_frequent_indels_pattern
4. check_type_of_mismatch


The option of the method
nuc_ref[str]:
read_info_ori[str]:
qualities[str]:
allele[str]:
allele_length[str]:
position[str]:
pileup_nb_reads[str]:
cut_off[str]:
minimum_coverage[str]:

return

"""
def pileup_extract_information(nuc_ref, read_info_ori, qualities, allele, allele_length, position, pileup_nb_reads,cut_off,minimum_coverage):
    ##=====================================================================================================##
    ## example: read_info_ori = "..,,,+1G...,-2GG,,^,^..^,,^A,^C,,,,$,$,"
    ## '.' => match read_forward, ',' => match read_reverse, '^' => start_of read, '$' => end of read, 'A'/'C'/'T'/'G' = > substitution, '*' deletion
    ## '^1'or other integers /'^A'/'^C'/'^G'/'^T' => mapping quality score
    ## '+1G' => insertion of G at next position, -2GG => deletion of two Gs at the next position
    ##=====================================================================================================##
    ## count '^' and '$' in read info ## high ratio => indicate large insertion or deletion if not located at begining/end of ref seq

    #nuc_ref -A, #read_info_ori -  ,,,,,,.,,,,,,,,...,,.,,,...,,,,..,,.,,,.,,,,,,,,.,,.,,,.,,,.^~.^~,^~,
    #qualities - FFFFIIFFIIIIIIIBFFIIFIIIFIIIIFIIFFIIFFFFFFBBFIIIFIFIFFIF<FFBBIB
    #allele - dfr_G, #allele_length - 498, #position - 1  # the posision of the read
    #pileup_nb_reads: 55 ,#cut_off: 84:50, #minimum_coverage: 5
    
    nb_start_end_info = len([elem for elem in list(read_info_ori) if re.match(r'[\$]',elem) != None or re.match(r'[\^]',elem) != None]) #count the number of ^ and $,  #nb_start_end_info  --> 55
    
    ratio_start_end = round(float(nb_start_end_info)*100/float(pileup_nb_reads),2) # percentage of ^ and $, compare to the total number of reads
    # ratio_start_end 100
    ## resub mapping quality score ^A,^T,^C,^G by '^', therefore keeping start of read info ##
    read_info = re.sub('[\^][\*ATCGatcg]','^',read_info_ori) # subtitute to '^' followed by ATCGatcg with ^, #read_info #  ^~,^~.^~,^~.^~,^
    ## uppercase string in reads_info and split into list (low case nucleotide may exist in reference)
    reads_info_list = list(read_info.upper()) #['^', '~', ',', '^', '~', '.',
    
    ## extract indel patterns [+/-][integer][ACTGxn] ##
    reads_info_modif_1 = combine_consecutive_integers_in_list(reads_info_list)
    #print "reads_info_modif_1", reads_info_modif_1 # reads_info_modif_1 ['^', '~', ',', '^', '~', '.
    read_parsed = combine_indels_pattern(reads_info_modif_1)
    #print "read_parsed", read_parsed  #read_parsed ['^', '~', ',', '^', '~', '.',-2AA,+1G,*,A,
    
    ## group read info per types (deletions,insertions,others)##
    deletions = [elem for elem in read_parsed if re.match(r'[-][\d]',elem) != None]
    #print "deletions", deletions # puts all the deletion in one array
    insertions = [elem for elem in read_parsed if re.match(r'[\+][\d]',elem) != None]
    #print "insertions", insertions # puts all the deletion in one array
    reads_per_types = [elem for elem in read_parsed if re.match(r'[\*\.,ACTGactg]',elem) != None] #* deletion in the previous base (*) and match or SNP
    
    ## filtering reads based on mapping quality and replacing reverse ',' by forward '.'##
    reads_high_quality = []
    for x,y in zip(reads_per_types,qualities): #NO QUALITY FOR DELETION, THEREFOR ONLY CALCULATE QUALITY FOR MATCH, SNP OR * ##qualities - FFFFIIFFIIIIIIIBFFIIFIIIFIIIIFIIFFIIFFFFFFBBFIIIFIFIFFIF<FFBBIB
        if ord(y) - 33 > 30:
            if x == ',': # if x is , replace it with .
                reads_high_quality.append('.')
            else:
                reads_high_quality.append(x)      #else keep  snp, or *
    nb_high_quality_reads = len(reads_high_quality)
    
    
    ## cut_off values by default:'84:50', can be modified using argument -c (cut_off)##
    ## first value cut off for match, substitutions => reported as not mixed if > 84 % indicate so
    ## second value cut off for indels => indels are reported if > 50% of reads indicate so
    values = cut_off.split(':')
    cut_off_top = round(float(values[0])/100,2)
    cut_off_bottom = 1.00 - cut_off_top
    indels_cut_off = round(float(values[1])/100,2)
    results = []
    mix_info = {}    
    
    if nb_high_quality_reads > 1:
        # count nb of indel patterns
        insertions_types = Counter(insertions) ## return dictionary=(key=type, value=nb_type), eg.{+1G=>54, +2GA=>1}
        ins_pattern_seq,ins_nb_pattern,ratio_pattern_insertion = extract_most_frequent_indels_pattern(insertions_types,pileup_nb_reads) # pileup_nb_reads = number of reads in the aligned
        deletion_types = Counter(deletions) ## return dictionary=(key=type,value=nb_type), eg.{-1C=>35}
        del_pattern_seq,del_nb_pattern,ratio_pattern_deletion = extract_most_frequent_indels_pattern(deletion_types,pileup_nb_reads)
        
        ## high quality reads are used at this stage ##
        
        ## count remaining nb of matchs, mismatchs, deletions(*) 
        total_reads_types = Counter(reads_high_quality) ## return dictionary=(key=type, value=nb_type), eg.{'.'=>75, 'A'=>2, 'C'=>1}
        nb_match = total_reads_types['.']
        mismatch_types = dict([(x,total_reads_types[x]) for x in total_reads_types if re.match(r'[ATCG]',x) != None])
        nb_mismatch = sum(mismatch_types.values())
        deletion_info = dict([(x,total_reads_types[x]) for x in total_reads_types if re.match(r'[\*]',x) != None])
        nb_deletion = sum(deletion_info.values())
        ## deduce ratio per type of info compared to total nb of reads at this position
        ratio_match = round(float(nb_match)/float(nb_high_quality_reads),2)
        ratio_mismatch = round(float(nb_mismatch)/float(nb_high_quality_reads),2)
        ratio_deletion = round(float(nb_deletion)/float(nb_high_quality_reads),2)
        
        ## return deduced info based on ratio/cut off value
        if ratio_match >= cut_off_top:
            results.extend(['match',nuc_ref, nb_match])
        elif ratio_match <= cut_off_bottom:
            if ratio_mismatch >= cut_off_top:
                predicted_nucleotide,nb_mismatch,extra_info = check_type_of_mismatch(mismatch_types,cut_off_top)
                mix_info = extra_info
                results.extend(['mismatch', predicted_nucleotide, nb_mismatch])
            elif ratio_deletion >= indels_cut_off:
                results.extend(['deletion','*',nb_deletion])
            elif ratio_mismatch < cut_off_top and ratio_deletion < indels_cut_off:
                results.extend(['mix','N',nb_high_quality_reads])                    
        else:
            if ratio_deletion > indels_cut_off:
                results.extend(['deletion','*',nb_deletion])
            else:
                if ratio_mismatch >= cut_off_bottom and ratio_deletion < indels_cut_off:
                    mismatch_nucleotide,nb_mismatch,extra_info = check_type_of_mismatch(mismatch_types,cut_off_top)
                    list_nucleotides = [nuc_ref,mismatch_nucleotide]
                    predicted_nuc = iupac_nucleotide_hash_code(list_nucleotides)
                    mix_info[nuc_ref] = ratio_match
                    mix_info[mismatch_nucleotide] = ratio_mismatch
                    results.extend(['mix',predicted_nuc, nb_high_quality_reads])
                else:
                    results.extend(['mix','N', nb_high_quality_reads])
        
        if ratio_pattern_insertion < indels_cut_off:
            results.extend(['ins', 0, 0])
        else:
            results.extend(['ins',ins_pattern_seq,ins_nb_pattern])
            
        if ratio_pattern_deletion < indels_cut_off:
            results.extend(['del', 0, 0])
        else:
            results.extend(['del',del_pattern_seq,del_nb_pattern])
        
        results.append(ratio_start_end)
        results.append(mix_info)  
    else:
        results.extend(['low_coverage','$',nb_high_quality_reads,'ins', 0, 0,'del', 0, 0])
        results.append(ratio_start_end)
        results.append(mix_info)
    ## return ['mapping_results', nuc, nb_reads,'ins',insert_seq,insert_nb,'del', del_seq, del_nb, ratio_start_end,mix_info]
    ## eg. ['mix','M',78,ins,0,0,del,0,0,0.1,{A=>0.62,C=>038}]
    ## eg. ['match','A',78,ins,0,0,del,0,0,0.1,{}]
    
    return results



"""
2.2
Called within the parsing_mpileup function

The option of the method[str]:
position_infos_ori[str]:
sequence_raw_ori[str]:
position_deletions[str]:
positions_with_accepted_depth[str]:
total_positions[str]:

return

"""
def update_infos_with_deletions(position_infos_ori,sequence_raw_ori,position_deletions,positions_with_accepted_depth,total_positions):
    ## update position with deletions eg.position_infos_ori = {12=>80,13=>70,14=>70}, sequence_raw_ori = {12=>A,13=>N,14=>N}, position_deletions = {12=>-2GG}
    ## update position_infos={12=>80,13=>*,14=>*}, sequence_raw={12=>A,13=>*,14=>*}
    positions_infos = position_infos_ori
    sequence_raw = sequence_raw_ori
    for position in position_deletions:
        if positions_with_accepted_depth.has_key(position):
            pattern_1 = position_deletions[position]
            nb_deleted_nuc = re.findall(r'[\d]{1,3}',pattern_1)
            range_start = int(position) + 1
            range_end = int(position)+ int(nb_deleted_nuc[0]) + 1
            nb_position_to_delete = range(range_start,range_end)
            for pos in nb_position_to_delete:
                sequence_raw[int(pos)] = '*'
                positions_infos[int(pos)] = '*'
    return positions_infos, sequence_raw 

"""
2.3
Called within the parsing_mpileup function

The option of the method:
dictionary[str]:
list_total_positions[str]:
positions_with_accepted_depth[str]:

return

"""
def update_infos_with_coverage(dictionary,list_total_positions,positions_with_accepted_depth):
    ## flag missing position eg. dictionary = {11=>A,13=>T,14=>C} to {11=>A,12=>"$",13=>T,14=>C}
    dictionary_updated = dictionary
    for tot_pos in list_total_positions: 
        if dictionary.has_key(tot_pos):
            if positions_with_accepted_depth.has_key(tot_pos):
                pass
            else:
                dictionary_updated[int(tot_pos)] = '$'
        else:
            dictionary_updated[int(tot_pos)] = '$'
    return dictionary_updated



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
2.4
Called within the parsing_mpileup function, uses pysam to pull missed bases
The option of the method:
fasta_file[str]:fasta_file /phengs/hpc_storage/home/rediat/python_code/gene_finder/gene_reference_data/staphylococcus-aureus-typing.1-0/reference.fasta
positions_with_accepted_depth_ori[str]:##positions_with_accepted_depth_ori {1: 55, 2: 55, 3: 56, 4:
position_mismatchs_ori[str]:#position_mismatchs_ori {}, #positions_indels_probabilities_ori {}
positions_indels_probabilities_ori[str]:# allele dfr_G, #allele_length 498, #sequence_raw_ori {1: 'A', 2: 'T', 3: 'G', 4: 'A', 5: 'A', 6:
allele[str]:
allele_length[str]:
sequence_raw_ori[str]:
sorted_bam[str]:#/home/rediat/gene_finder_fastq/29188_H14276008707-939/GF_output/tmp/29188_H14276008707-939-1.sorted.bam

return
"""
def complete_missing_seq(fasta_file,positions_with_accepted_depth_ori,position_mismatchs_ori,positions_indels_probabilities_ori,allele,allele_length,sequence_raw_ori,sorted_bam):
    
    positions_with_accepted_depth = positions_with_accepted_depth_ori #????
    position_mismatchs = position_mismatchs_ori #????
    positions_indels_probabilities = positions_indels_probabilities_ori #????
    sequence_raw = sequence_raw_ori #????
    target_biosequence_list = []           
    for key in sorted(sequence_raw.iterkeys()):
        target_biosequence_list.append(sequence_raw[key])   
    
    last_mapped_position =  max(sorted(positions_with_accepted_depth.keys())) #{1: 55, 2: 55, 3: 56, 4:, last_mapped_position # 498
    first_mapped_position =  min(sorted(positions_with_accepted_depth.keys())) #first_mapped_position # 1

    if first_mapped_position != 1 and first_mapped_position < 20 and float(len(positions_with_accepted_depth.keys())) > float(allele_length*0.8):
        target_biosequence_string = ''.join([elem for elem in target_biosequence_list])
        pattern_search = re.findall(r'[ACTG]{7,10}',target_biosequence_string)[0] #????#first 7 to 10 consecutive nuc in predicted sequence to search reads#
        pattern_search_index = target_biosequence_string.find(pattern_search)
        missing_positions = pattern_search_index 
        if positions_indels_probabilities.has_key(first_mapped_position):
            del positions_indels_probabilities[first_mapped_position]
        missing_nuc, nb_reads = pysam_search(sorted_bam,allele,pattern_search,pattern_search_index,pattern_search_index + 10, missing_positions)
        if missing_nuc != 'nil':
            add_mism = check_missing_nuc_mismatch(fasta_file,allele,0,pattern_search_index +1,missing_nuc)
            added_positions = range(1,missing_positions+1)
            add_mism_checked = dict([(x,add_mism[x]) for x in add_mism if add_mism[x][0] != add_mism[x][1]])
            position_mismatchs.update(add_mism_checked)
            for x in added_positions:
                positions_with_accepted_depth[x] = nb_reads
                sequence_raw[x] = list(missing_nuc)[x-1]
        else:
            pass
    ## add sequences at end of gene if missing  
    if last_mapped_position != allele_length and last_mapped_position > allele_length - 20 and float(len(positions_with_accepted_depth.keys())) > float(allele_length*0.8):
        target_biosequence_string = ''.join([elem for elem in target_biosequence_list])
        pattern_search = re.findall(r'[ACTG]{7,10}',target_biosequence_string)[-1] #last 7 to 10 consecutive nuc in predicted sequence to search reads#
        pattern_search_index_start = [x.start() for x in re.finditer(pattern_search,target_biosequence_string)][-1]
        pattern_search_index_end = pattern_search_index_start + len(pattern_search)
        missing_positions = allele_length - pattern_search_index_end
        if positions_indels_probabilities.has_key(last_mapped_position):
            del positions_indels_probabilities[last_mapped_position]
        missing_nuc, nb_reads = pysam_search(sorted_bam,allele,pattern_search,pattern_search_index_end,pattern_search_index_end -10,missing_positions)
        if missing_nuc != 'nil':
            add_mism = check_missing_nuc_mismatch(fasta_file,allele,pattern_search_index_end,allele_length,missing_nuc)
            add_mism_checked = dict([(x,add_mism[x]) for x in add_mism if add_mism[x][0] != add_mism[x][1]])
            position_mismatchs.update(add_mism_checked)                            
            added_positions = range(pattern_search_index_end+1,pattern_search_index_end+1+missing_positions)
            combined_positions_missing_nuc = zip(added_positions,list(missing_nuc))
            for y in combined_positions_missing_nuc:
                sequence_raw[y[0]] = y[1]
            for x in added_positions:
                positions_with_accepted_depth[x] = nb_reads
        else:
            pass
    return positions_with_accepted_depth,position_mismatchs,sequence_raw,positions_indels_probabilities

"""
#2.4.1 called within the complete_missing_seq which is called by parsing_mpileup function

The option of the method:
sam_file[str]:
allele[str]:
pattern[str]:
map_1[str]:
map_2[str]:
missing_positions[str]:

return
"""
def pysam_search(sam_file,allele,pattern,map_1,map_2,missing_positions):
    length_pattern = len(pattern)
    seq_to_add = []
    seq_mapped_read = []
    samfile = pysam.Samfile(sam_file,'rb')
    if map_1 < map_2:
        for alignedread in samfile.fetch(allele,map_1,map_2):
            read_seq = alignedread.seq
            seq_mapped_read.append(read_seq)
            samfile.close()
    else:
        for alignedread in samfile.fetch(allele,map_2,map_1):
            read_seq = alignedread.seq
            seq_mapped_read.append(read_seq)
            samfile.close()
    for seq in seq_mapped_read:
        match_postion = seq.find(pattern)
        if match_postion == -1:
            pass
        else:
            if map_1 < map_2 :
                if match_postion - missing_positions >= 0:
                    substring = seq[match_postion - missing_positions :match_postion]
                    seq_to_add.append(str(substring))
                else:
                    substring = seq[0:match_postion]
                    seq_to_add.append(str(substring))
            else:
                substring = seq[match_postion + length_pattern :match_postion + length_pattern + missing_positions]
                seq_to_add.append(str(substring))  
    missing_seq = Counter(seq_to_add)
    if len(missing_seq) != 0:
        nuc_seq = (missing_seq).most_common(1)[0][0]
        if len(nuc_seq) == missing_positions:
            ## nb of reads from which missing sequences were extracted
            nb_reads = missing_seq[nuc_seq]
            return nuc_seq, nb_reads
        else:
            nb_missing = missing_positions - len(nuc_seq)
            extension = ['N'] * nb_missing
            extracted = list(nuc_seq)
            extracted.extend(extension)
            nb_reads = missing_seq[nuc_seq]
            new_nuc_seq = ''.join(map(str,extracted))
            return new_nuc_seq, nb_reads
    else:
        return 'nil', 'nil'
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.5.called with parsing_mpileup
def updated_infos_with_insertions(positions_with_accepted_depth,position_insertions,sequence_raw_ori):
    ## update position with insertions eg.sequence_raw_ori = {12=>A,13=>T,14=>G}, position_deletions = {12=>+2AT}
    ## update sequence_raw={12=>AAT,13=>T,14=>G}
    total_inserted_nuc = []
    sequence_raw = sequence_raw_ori
    for position in position_insertions:
        if positions_with_accepted_depth.has_key(position):
            pattern_2 = position_insertions[position]
            insertion = re.findall(r'[ACGT]+',pattern_2)
            new_value = sequence_raw[position] + insertion[0]
            sequence_raw[position] = str(new_value)
            total_inserted_nuc.append(len(insertion[0]))
    return sequence_raw, total_inserted_nuc



###########################################################################################################################



#2.1.1.called within pileup_extract_information function
def combine_consecutive_integers_in_list(list_of_elem):#list_of_elem: reads_info_modif_1 ['^', '~', ',', '^', '~', '.	
    ## Replace two consecutive integer elements into one element in a list (eg. [a,a,c,1,2,d,e,1,f] => [a,a,c,12,d,e,1,f]
    # processing indel there are indeles in the read the will be split in the array, i.e 12 will be split 12acct
    output = []
    i = 0
    while len(list_of_elem) > 1:
        a = list_of_elem[i]
        b = list_of_elem[i+1]
        if is_integer(a) and is_integer(b):
            nbs = ''.join(list_of_elem[i:i+2])
            list_of_elem[i:i+2] = []
            list_of_elem.insert(0,nbs)
        else:
            output.append(a)
            del list_of_elem[i]
    output.extend(list_of_elem)
    
    return output

#2.1.2.called within pileup_extract_information function
def combine_indels_pattern(list_of_elem):
    ## indel_pattern = [+/-][interger][integer x letter] (eg.-2GG or +1C)
    ## to be used after method:combine_consecutive_integers_in_list
    ## combine indel pattern in the list (eg. [a,a,+,1,c,a,a,-,10,c,c,c,c,c,c,c,c,c,c,a,a,a] => [a,a,+1c,a,a,-10cccccccccc,a,a,a]
    index_numbers = []
    for (index,item) in enumerate(list_of_elem):
        if is_integer(item) == True:
            index_numbers.append(index)
    for i in reversed(index_numbers):
        if list_of_elem[i-1] == '+' or list_of_elem[i-1] == '-':
            extension = int(list_of_elem[i]) + 1
            new_element = ''.join(list_of_elem[i-1:i+extension])
            list_of_elem[i-1 :i+extension] = [str(new_element)]
    return list_of_elem



#2.1.3.called within pileup_extract_information function
def extract_most_frequent_indels_pattern(indels_types,pileup_nb_reads): #?? why the most frequent pattern nedded????
    ## return a list of three value (pattern,count,ratio_insertion) (eg. +1G,50,0.98)
    if len(indels_types) == 0:
        return 0,0,0
    else:
        for pattern,count in indels_types.most_common(1):      
            ratio_insertion = round(float(int(count)/int(pileup_nb_reads)),2)
            return pattern,count,ratio_insertion
	
#----------------------------------------------------------------------
#2.1.4.called within pileup_extract_information function
def check_type_of_mismatch(nb_reads_types,cut_off):
    ## input is a dictionary (key=mismatch_type, value=nb_mismatch_type), eg. {A=>24, C=>28)
    extra_info = {}
    totat_nb_mismatch = sum(nb_reads_types.values())
    mis_nuc = sorted(nb_reads_types,key=nb_reads_types.get, reverse=True)
    if nb_reads_types[mis_nuc[0]] >= (totat_nb_mismatch * cut_off):
        return mis_nuc[0],nb_reads_types[mis_nuc[0]],extra_info
    else:
        mix_mismatch_nb = nb_reads_types[mis_nuc[0]] + nb_reads_types[mis_nuc[1]]
        if mix_mismatch_nb >= (totat_nb_mismatch * cut_off):
            ## take the two_most_commun_mismatchs = [mis_nuc[0],mis_nuc[1]] and return iupac code eg.( amix of A and C return M)
            predicted_nuc = iupac_nucleotide_hash_code([mis_nuc[0],mis_nuc[1]])
            ratio_0 = round(float(nb_reads_types[mis_nuc[0]]) / float(mix_mismatch_nb),2)
            ratio_1 = round(float(nb_reads_types[mis_nuc[1]]) / float(mix_mismatch_nb),2)
            extra_info[mis_nuc[0]] = ratio_0
            extra_info[mis_nuc[1]] = ratio_1
            return predicted_nuc,mix_mismatch_nb,extra_info
        else:
            return 'N', totat_nb_mismatch,extra_info


# called within the 2.1.4 check_type_of_mismatch function
"""
if the positio shows a mixed position, iupac bases are assigned. e.g. If the mixed base are (A and G), sub with R


"""
def iupac_nucleotide_hash_code(list_two_nucleotide):
	
    code = []
    list_two_nucleotide.sort()
    if list_two_nucleotide[0] == 'A':
        if list_two_nucleotide[1] == 'G':
            code.append('R')
        elif list_two_nucleotide[1] == 'T':
            code.append('W')
        elif list_two_nucleotide[1] == 'C':
            code.append('M')
        else:
            code.append('N')
    elif list_two_nucleotide[0] == 'C':
        if list_two_nucleotide[1] == 'T':
            code.append('Y')
        if list_two_nucleotide[1] == 'G':
            code.append('S')
        else:
            code.append('N')
    elif list_two_nucleotide[0] == 'G':
        if list_two_nucleotide[1] == 'T':
            code.append('K')
        else:
            code.append('N')
    else:
        code.append('N')
    return code[0]


#------------------------------------------------------------------------------
#called by 2.4 complete_missing_seq
def check_missing_nuc_mismatch(fasta_file,allele,first_position,last_position,missing_nuc):
    report = {}
    comp_list = []
    for ref in SeqIO.parse(fasta_file,'fasta'):
        if ref.id == allele:
            reference_seq_list = list(ref.seq)
            if first_position == 0:
                comp_list = zip(reference_seq_list[first_position:last_position-1],list(missing_nuc))
            else:
                comp_list = zip(reference_seq_list[first_position:last_position],list(missing_nuc))
            for (index,element) in enumerate(comp_list):
                if first_position == 0:
                    report[index + 1] = [element[0],element[1]]
                else:
                    report[first_position + index +1] = [element[0],element[1]]
    return report







