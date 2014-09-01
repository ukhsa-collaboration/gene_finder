#!/usr/bin/env python

from __future__ import division
import logging, pprint
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
import shutil

module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
	module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
	if module_folder not in sys.path:
		sys.path.insert(1, module_folder)

import gene_finder_functions
import generate_mpileup_file
import parsing_mpileup
import extract_quality_metrics
import log_writer
from utility_functions import *




"""
Function
1.Index reference file 
2.Generating SAM, BAMs and mpileup
3.Parse a pileup file 
4.extract_quality_metrics
5.Find best_hit
6.report result in xml file

The option of the method
outdir[str]: path to output_directory
fasta_file[str]: full path to reference fasta file
fastq_files[list]: full path to the fastq file
bowtie_options[list]: list bowtie_options, default:['-q', '--very-sensitive-local', '--no-unal', '-a']
ids[str]: sample id (eg. NGSLIMSID_molis id)
cut_off[str]: first number is the cut off to be used used to identify mix at each position when parsing the pileup file (e.g. ......,,,,,,A,T,,,,, if match >= 84 => no mix else mix)
	second number is the cut off to be used used to identify deletions at each position when parsing the pileup file (e.g. ......,,,*****,*,,,,, if nb of * >= 50 => deletion)
minimum_coverage[str]: check if the min  depth of coverage is 5% of the maximum depth of coverage 
log_directory[str]: full path to the log file
workflow_name[str]: species_workflow (eg. staphylococcus_typing)
version[str]: verion number (eg. 1-0-0)

Returns
Write result in xml file
"""

def run_gene_finder(outdir,fasta_file,workflow_file, fastq_files,bowtie_options,ids,cut_off,minimum_coverage,log_directory, workflow_name="",version=""):

    stderr_log_output = log_directory + "/" + 'gene_finder'+ ".stderr"
    stdout_log_output = log_directory + "/" + 'gene_finder'+ ".stdout"
    logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)
    
    forward_fastq = fastq_files[0]
    reverse_fastq = fastq_files[1]
    if not os.path.exists(outdir + '/tmp'):
        os.mkdir(outdir + '/tmp')
	
    path_to_tmp_file = outdir + '/tmp'    
    shutil.copy(fasta_file, path_to_tmp_file)
    fasta_file = path_to_tmp_file + "/reference.fasta"
    
    try_and_except(stderr_log_output,generate_mpileup_file.index_reference,fasta_file,logger)
    try_and_except(stderr_log_output,generate_mpileup_file.samtools_faidx,fasta_file,logger)
    try_and_except(stderr_log_output, generate_mpileup_file.generate_mpileup,path_to_tmp_file,fasta_file,forward_fastq,reverse_fastq,outdir,workflow_name,version,ids,bowtie_options,logger)
    pileup_dictionary = try_and_except(stderr_log_output, parsing_mpileup.read_mipelup,fasta_file,cut_off,minimum_coverage,outdir,workflow_name,version,ids,logger)
    pileup_quality_dictionary= try_and_except(stderr_log_output,extract_quality_metrics.extract_quality_metrics, pileup_dictionary)
    path_to_reference_folder = os.path.dirname(fasta_file)
    xml_output,control_coverage_value= try_and_except(stderr_log_output,best_hit,workflow_file,pileup_quality_dictionary,path_to_reference_folder,outdir)
    try_and_except(stderr_log_output,report_xml,outdir,workflow_name,version,ids,xml_output,control_coverage_value)
    



"""
Function
If the length of target and reference sequence are equal:
1. translate the nuclotide sequence
2. compare the protein sequence between the reference  and target sequence

The option of the method   
reference_seq[<class 'Bio.Seq.Seq'>]:reference nuclotide  sequence .
target_fasta_sequence[<class 'Bio.Seq.Seq'> ]:target nuclotide  sequence

Returns
mutations[dict]: Position and change of  mutated amino acid sequence. Also return region (whether the change is in non stop or early stop codon region)
e.g.  {145: 'no_stop', 37: 'I-X'}
"""
def compare_protein_sequence(reference_seq,target_seq):
	
    mutations = {}
    ref_prot = list(reference_seq.translate())
    target_prot = list(target_seq.translate())
    comparative = zip(ref_prot,target_prot)
    for index,item in enumerate(comparative):
        if item[0] == item[1]:
            pass
        else:
            report = ''.join(map(str,[item[0],'-',item[1]]))
            mutations[(index+1)] = report
    try:
        stop_codon_pos = list(target_prot).index('*')
        if stop_codon_pos != (len(list(target_prot)) - 1):
            mutations[stop_codon_pos +1] = 'early_stop'
    except ValueError:
        mutations[len(list(target_prot))] = 'no_stop'
    return mutations

"""
Function
If the length of target and reference sequence are not equal:
1. translate the nuclotide sequence
2. Compare the protein sequence between the reference  and target sequence using Clustalw. 

The option of the method   
reference_seq[<class 'Bio.Seq.Seq'>]:  reference nuclotide  sequence .
target_fasta_sequence[<class 'Bio.Seq.Seq'> ]:  target nuclotide  sequence
outdir[str]: full path to the output dir  

Returns
final_report[dict]: Position and change of  mutated amino acid sequence. Also return region (whether the change is in non stop or early stop codon region)
e.g.  {145: 'no_stop', 37: 'I-X'}
"""
   
def clustalw_align(reference_seq,target_seq,output):

    clustal_output = output + '/tmp/clustal_align'
    if os.path.exists(clustal_output):
        pass
    else:
        os.mkdir(clustal_output)
    reference_prot = reference_seq.translate()
    target_prot = target_seq.translate()
    multi_fasta = '>reference' + '\n' + reference_prot + '\n' +'>target' + '\n' + target_prot + '\n'
  
    with open(clustal_output + '/multi_fasta.fasta','w') as my_file:
        my_file.write(str(multi_fasta))
        my_file.close()
    clustalw_cline = ClustalwCommandline('clustalw2', infile= clustal_output + '/multi_fasta.fasta')
    stdout, stderr = clustalw_cline()
    alignment = AlignIO.read(clustal_output + '/multi_fasta.aln', "clustal")
    reference_alig_seq,target_align_seq = [record.seq for record in alignment]
    final_report = {}
    counter = 0
    for index,item in enumerate(zip(reference_alig_seq,target_align_seq)):
        if item[0] == item[1]:
            pass
        else:
            if item[0] != item[1]:
                if item[0] == '-':
                    report = ''.join(map(str,[item[0],'-',item[1]]))
                    final_report[(index-counter)] = report
                    counter += 1
                else:
                    report = ''.join(map(str,[item[0],'-',item[1]]))
                    final_report[(index+1-counter)] = report
    try:
        stop_codon_pos = list(target_prot).index('*')
        if stop_codon_pos != (len(list(target_prot)) - 1):
            final_report[stop_codon_pos +1] = 'early_stop'
    except ValueError:
        final_report[len(list(target_prot))] = 'no_stop'
    return final_report

"""
Function
Compare the protien sequence between the reference  and target sequence. 

The option of the method   
reference_seq[<class 'Bio.Seq.Seq'>]: reference nuclotide  sequence .
target_fasta_sequence[<class 'Bio.Seq.Seq'> ]: target nuclotide  sequence
outdir[str]: full path to the output dir  

Returns
mutations[dict]Position and change of  mutated amino acid sequence. Also return region (whether the change is in non stop or early stop codon region)
e.g.  {145: 'no_stop', 37: 'I-X'}
"""
   
def check_variant_transcription(reference_seq,target_fasta_sequence,outdir):
  
    if len(target_fasta_sequence) == 1:
        if len(target_fasta_sequence[0].seq) == len(reference_seq):           
            mutations = compare_protein_sequence(reference_seq,target_fasta_sequence[0].seq)
            return mutations
        else:
            target_ref_sequence = target_fasta_sequence[0].seq
            mutations = clustalw_align(reference_seq,target_ref_sequence,outdir)
            return mutations
    else:
        return "ND"

"""
Function

The option of the method   
pileup_hash[dict]: primary  key = allele names, secondary key= type of quality metrics  (e.g.:housekeeping_detection,ref_cut_off, coverage_distribution,position_deletions,protein_alterations,coverage etc}   andvalue = value assigned to the secondary key
path_to_reference_folder[str]: full path to the refrence folder 
outdir[str]: full path to the output dir 

Returns
xml_output[dict]:
primary  key = allele names, secondary key= type of quality metrics   (e.g.:housekeeping_detection,ref_cut_off, coverage_distribution,position_deletions,protein_alterations,coverage etc}   and
value = value assigned to the secondary key
	  
control_coverage_value[str]:

""" 
def best_hit(workflow_file,pileup_hash,path_to_reference_folder,outdir):
	
    control_coverage = []
    list_of_keys = pileup_hash.keys()
    xml_output = {}
    with open(workflow_file) as workflow_info:
    #with open(path_to_reference_folder + '/workflow.txt') as workflow_info:
        for line in workflow_info:
            final_output = []
            gene,ref_homology,workflow,start_pos,end_pos = line.split()
            if any(key.split('_')[0] == gene for key in pileup_hash):
                best_hits = []
                gene_by_group = [elem for elem in list_of_keys if elem.split('_')[0] == gene]
                results_by_coverage = {}
                results_by_mismatch = {}
                results_by_insertion = {}
                results_by_mix = {}
                for variant in gene_by_group:
                    key_coverage = pileup_hash[variant]['coverage']
                    key_mismatch = len(pileup_hash[variant]['position_nuc_mismatchs'])
                    key_insertion = len(pileup_hash[variant]['position_insertions'])
                    if results_by_insertion.has_key(key_insertion):
                        results_by_insertion[key_insertion].append(variant)
                    else:
                        results_by_insertion[key_insertion] = [variant]
                    if results_by_coverage.has_key(key_coverage):
                        results_by_coverage[key_coverage].append(variant)
                    else:
                        results_by_coverage[key_coverage] = [variant]
                    if results_by_mismatch.has_key(key_mismatch):
                        results_by_mismatch[key_mismatch].append(variant)
                    else:
                        results_by_mismatch[key_mismatch] = [variant]                        
                sorted_max_coverage = sorted(results_by_coverage, key = lambda x:float(x), reverse=True)
                sorted_mismatch_results = sorted(results_by_mismatch)
                sorted_results_insertion = sorted(results_by_insertion)
                
                for max_coverage in sorted_max_coverage:
                    genes_max_coverage = results_by_coverage[max_coverage]
                    for mismatch_nb in sorted_mismatch_results:
                        genes_min_mismatch = results_by_mismatch[mismatch_nb]
                        common_genes = [elem for elem in genes_min_mismatch if elem in genes_max_coverage]
                        if len(common_genes) != 0:
                            for nb_insertion in results_by_insertion:
                                genes_min_insertion = sorted(results_by_insertion[nb_insertion])
                                best_hit = [elem for elem in common_genes if elem in genes_min_insertion]
                                if len(best_hit) != 0:
                                    best_hits.append(best_hit[0])
                                    break
                                else:
                                    pass
                        else:
                            pass
                        
                target_fasta_sequence = pileup_hash[best_hits[0]]['predicted_seq']
                hash_reference_seq = {}
                for ref in SeqIO.parse(path_to_reference_folder + '/reference.fasta','fasta'):
                    hash_reference_seq[ref.id] = ref.seq                   
                if workflow == 'mutant':
                    control_coverage.append(float(pileup_hash[best_hits[0]]['depth'].split(':')[0]))
                    mutation_info = {}
                    reference_seq = hash_reference_seq[best_hits[0]]
                    if pileup_hash[best_hits[0]]['coverage'] == 100.0:
                        if len(target_fasta_sequence) == 1:
                            if len(target_fasta_sequence[0].seq) == len(reference_seq):
                                mutations = compare_protein_sequence(reference_seq,target_fasta_sequence[0].seq)
                                mutation_info['positions_prot_modifications'] = mutations
                            else:
                                target_ref_sequence = target_fasta_sequence[0].seq
                                mutations = clustalw_align(reference_seq,target_ref_sequence,outdir)
                                mutation_info['positions_prot_modifications'] = mutations
                        else:
                            mutation_info['positions_prot_modifications'] = "ND"
                    else:
                        mutation_info['positions_prot_modifications'] = "ND"
                        
                    if mutation_info['positions_prot_modifications'] != "ND":
                        if len(mutation_info['positions_prot_modifications']) != 0:
                            relevent_position = dict([(key,mutation_info['positions_prot_modifications'][key]) for key in mutation_info['positions_prot_modifications'] if int(key) >= int(start_pos) and int(key) <= int(end_pos)])
                            pileup_hash[best_hits[0]]['relevent_position'] = relevent_position
                        else:
                            pileup_hash[best_hits[0]]['relevent_position'] = {}
                    else:
                        pileup_hash[best_hits[0]]['relevent_position'] = "ND"
                        
                    pileup_hash[best_hits[0]]['ref_cut_off'] = int(ref_homology)
                    xml_output[best_hits[0]] = pileup_hash[best_hits[0]]
                elif workflow == 'variant':
                    reference_seq = hash_reference_seq[best_hits[0]]
                    if pileup_hash[best_hits[0]]['coverage'] == 100: #and best_hits[0].split("_")[0] != "23s":
                        alterations = check_variant_transcription(reference_seq,target_fasta_sequence,outdir)
                        pileup_hash[best_hits[0]]['protein_alterations'] = alterations
                    else:
                        pileup_hash[best_hits[0]]['protein_alterations'] = 'ND'
                    pileup_hash[best_hits[0]]['ref_cut_off'] = int(ref_homology)
                    xml_output[best_hits[0]] = pileup_hash[best_hits[0]]
                elif workflow == 'regulator':
                    pileup_hash[best_hits[0]]['nucleotide_alterations'] = {}
                    total_position_mismatch = pileup_hash[best_hits[0]]['position_nuc_mismatchs']
                    continous_coverage = []
                    for section_covered in pileup_hash[best_hits[0]]['coverage_distribution']:
                        coverage_details = section_covered.split('-')
                        if int(start_pos) >= int(coverage_details[0]) and int(end_pos) <= int(coverage_details[1]):
                            continous_coverage.append("true")
                    if len(continous_coverage) == 1:
                        if len(total_position_mismatch) != 0:
                            for sp_position in total_position_mismatch:
                                if int(sp_position) >= int(start_pos) and int(sp_position) <= int(end_pos):
                                    pileup_hash[best_hits[0]]['nucleotide_alterations'][sp_position] = total_position_mismatch[sp_position]
                    else:
                        pileup_hash[best_hits[0]]['nucleotide_alterations'] = 'ND'
                    pileup_hash[best_hits[0]]['ref_cut_off'] = int(ref_homology)
                    xml_output[best_hits[0]] = pileup_hash[best_hits[0]]
            else:
                if workflow == 'mutant':
                    xml_output[gene] = {}
                    xml_output[gene]['housekeeping_detection'] = 'ND'
                else:
                    pass           
    if len(control_coverage) > 0:       
        control_coverage_value = sorted(control_coverage)[0]
    else:
        control_coverage_value = 'failed'
    
    return xml_output,control_coverage_value


"""
Function
Report result in xml file 

The option of the method   
outdir[str]: full path to the output dir
workflow_name[str]: species_workflow (eg. staphylococcus_typing
version[str]: version number
prefix[str]: sampleid (e.g. NGSLIMS_molis number, 9188_H14276008707)
xml_output[dict]:primary key=allele name, secondary key = type of quality metrics and value= value of quality metrics
control_coverage_value[str]:

Returns
print result in xml file

""" 
def report_xml(outdir,workflow_name,version,prefix,xml_output,control_coverage_value):

    genome_id = prefix.split(".")[0]   
    xml_log_file = open(outdir + "/" + genome_id + ".results.xml", "w")
    root = etree.Element("ngs_sample", id = genome_id)
    workflow = etree.SubElement(root, "workflow", value=workflow_name, version = version)
    coverage_control = etree.SubElement(root, "coverage_control", value=str(control_coverage_value))
    results = etree.SubElement(root, 'results')
    
    for gene_cluster, value in xml_output.items():
        #gene = gene_cluster.split('_')[0]
	gene = gene_cluster
        result = etree.SubElement(results, "result", type="gene", value = gene)
        if value.has_key('housekeeping_detection'):
            etree.SubElement(result, "result_data", type="mode", value = 'mutant')
            etree.SubElement(result, "result_data", type="detection", value = 'not_detected')
        else:        
            if value.has_key('relevent_position'):
                etree.SubElement(result, "result_data", type="mode", value = 'mutant')
                if value['relevent_position'] != "ND":
                    if len(value['relevent_position']) != 0:
                        report = []
                        for position,mutation in sorted(value['relevent_position'].items()):
                            relevent_info =  str(position) + ':' + mutation
                            report.append(relevent_info)
                        mix = ';'.join(report)
                        etree.SubElement(result, "result_data", type="mutation", value = mix)
                    else:
                        etree.SubElement(result, "result_data", type="mutation", value = 'None')
                else:
                    etree.SubElement(result, "result_data", type="mutation", value = 'ND')
            elif value.has_key('nucleotide_alterations'):
                etree.SubElement(result, "result_data", type="mode", value = 'regulator')
                #print value['nucleotide_alterations']
                if value['nucleotide_alterations'] != "ND":
                    if len(value['nucleotide_alterations']) != 0:
                        report = []
                        for position,mutation in sorted(value['nucleotide_alterations'].items()):
                            relevent_info =  str(position) + ':' + '-'.join(mutation)
                            report.append(relevent_info)
                        mix = ';'.join(report)
                        etree.SubElement(result, "result_data", type="modifications", value = mix)
                    else:
                        etree.SubElement(result, "result_data", type="modifications", value = 'None')
                else:
                    etree.SubElement(result, "result_data", type="modifications", value = 'ND')                 
            else:
                etree.SubElement(result, "result_data", type="mode", value = 'variant')
                if value['protein_alterations'] != "ND":
                    if len(value['protein_alterations']) != 0:
                        report = []
                        for position,mutation in sorted(value['protein_alterations'].items()):
                            relevent_info =  str(position) + ':' + mutation
                            report.append(relevent_info)
                        mix = ';'.join(report)
                        etree.SubElement(result, "result_data", type="alterations", value = mix)
                    else:
                        etree.SubElement(result, "result_data", type="alterations", value = 'None')
                else:
                    etree.SubElement(result, "result_data", type="alterations", value = 'ND')
            
            if value['homology'] >= float(value['ref_cut_off']) and value['coverage'] == float(100):
                etree.SubElement(result, "result_data", type="detection", value = 'present')
            elif value['homology'] >= float((value['ref_cut_off'])) and value['coverage'] != float(100) or value['homology'] < float((value['ref_cut_off'])) and value['coverage'] == float(100):
                etree.SubElement(result, "result_data", type="detection", value = 'uncertain')
                all_contigs = []
                #print value['predicted_seq']
                for contig_predicted in value['predicted_seq']:
                    all_contigs.append(str(contig_predicted.id))
                    all_contigs.append(str(contig_predicted.seq))
                fasta_contigs_file = open(outdir + '/' + gene + '_contigs.fasta',"w")
                fasta_contigs_file.write('\n'.join(all_contigs))
                fasta_contigs_file.close 
            else:
                etree.SubElement(result, "result_data", type="detection", value = 'absent')
            
            etree.SubElement(result, "result_data", type="Coverage", value = str(value['coverage']))
            etree.SubElement(result, "result_data", type="Homology", value = str(value['homology']))
            etree.SubElement(result, "result_data", type="Depth", value = str(value['depth']))
            etree.SubElement(result, "result_data", type="coverage_distribution", value = str(value['coverage_distribution']))
            if len(value['positions_mix']) > 0:
                report = []
                for position,mutation in sorted(value['positions_mix'].items()):
                    all_values = mutation.values()
                    new = [str(elem) for elem in all_values]
                    mismatch_info = str(position) + ':' + '-'.join(mutation.keys()) + '-' + '-'.join(new)
                    report.append(mismatch_info)
                mix = ';'.join(report)
                etree.SubElement(result, "result_data", type="Mix", value = mix)
            else:
                etree.SubElement(result, "result_data", type="Mix", value = 'None')
            
            if len(value['probability_big_indels']) > 0:
                etree.SubElement(result, "result_data", type="large_indels", value = str(sorted(value['probability_big_indels'].items())))
            else:
                etree.SubElement(result, "result_data", type="large_indels", value = 'None')
            
            if len(value['position_nuc_mismatchs']) > 0:
                report = []
                for position,mutation in sorted(value['position_nuc_mismatchs'].items()):
                    mismatch_info = str(position) + ':' + '-'.join(map(str,mutation))
                    report.append(mismatch_info)
                mix = ';'.join(report)
                etree.SubElement(result, "result_data", type="mismatch", value = mix)
            else:
                etree.SubElement(result, "result_data", type="mismatch", value = 'None')
            #etree.SubElement(result, "result_data", type="sequence", value = str(value['predicted_seq']))
    print etree.tostring(root, pretty_print=True)
    print >> xml_log_file, etree.tostring(root, pretty_print=True)
   


    
