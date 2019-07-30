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


def run_gene_finder(logger,outdir,fasta_file,workflow_file, fastq_files,bowtie_options,ids,cut_off,minimum_coverage,log_directory, workflow_name="",version=""):
    """
    Function: run gene finder:
    - create tmp file in gene_finder results folder
    - index reference file 
    - generate SAM, BAMs and mpileup
    - Parse a pileup file 
    - extract_quality_metrics
    - find best_hits
    - report result in xml file
    
    Input:
    - outdir[str]: path to output_directory
    - fasta_file[str]: full path to reference fasta file
    - fastq_files[list of two elements]: full path to the fastq_forward and fastq_reverse
    - bowtie_options[list]: list bowtie_options, default:['-q', '--very-sensitive-local', '--no-unal', '-a']
    - ids[str]: sample id (eg. NGSLIMSID_molis id)
    - cut_off[str]: first number is the cut off to be used used to identify mix at each position when parsing the pileup file (e.g. ......,,,,,,A,T,,,,, if match >= 84 => no mix else mix)
        second number is the cut off to be used used to identify deletions at each position when parsing the pileup file (e.g. ......,,,*****,*,,,,, if nb of * >= 50 => deletion)
    - minimum_coverage[integer]: default [5] coverage at any position should be over 5% of the maximum coverage depth for the gene to be considered acceptable otherwise gap.
    - log_directory[str]: full path to the logs folder
    - workflow_name[str]: species_workflow (eg. staphylococcus_typing)
    - version[str]: verion number (eg. 1-0-0)
    
    Returns
    Write result in xml file
    """
    stderr_log_output = log_directory + "/" + 'gene_finder'+ ".stderr"    
    forward_fastq = fastq_files[0]
    reverse_fastq = fastq_files[1]
    if not os.path.exists(outdir + '/tmp'):
        os.mkdir(outdir + '/tmp')
    
    path_to_tmp_file = outdir + '/tmp'    
    shutil.copy(fasta_file, path_to_tmp_file)
    fasta_file = path_to_tmp_file + "/reference.fasta"
    
    try_and_except(stderr_log_output,generate_mpileup_file.index_reference,fasta_file,logger)
    try_and_except(stderr_log_output,generate_mpileup_file.samtools_faidx,fasta_file,logger)
    try_and_except(stderr_log_output,generate_mpileup_file.generate_mpileup,path_to_tmp_file,fasta_file,forward_fastq,reverse_fastq,outdir,workflow_name,version,ids,bowtie_options,logger)
    pileup_dictionary = try_and_except(stderr_log_output, parsing_mpileup.read_mipelup,fasta_file,cut_off,minimum_coverage,outdir,workflow_name,version,ids,logger)
    pileup_quality_dictionary= try_and_except(stderr_log_output,extract_quality_metrics.extract_quality_metrics, pileup_dictionary)
    path_to_reference_folder = os.path.dirname(fasta_file)
    xml_output,control_coverage_value,mix_indicator= try_and_except(stderr_log_output,best_hit,workflow_file,pileup_quality_dictionary,path_to_reference_folder,outdir)
    try_and_except(stderr_log_output,report_xml,outdir,workflow_name,version,ids,xml_output,control_coverage_value,mix_indicator)
    for root, dirs, files in os.walk(path_to_tmp_file):
        print "dirs", dirs
        for currentFile in files:
            exts=('.sam', '.mod', '.bt2','fai')
            if any(currentFile.lower().endswith(ext) for ext in exts):
                os.remove(os.path.join(root, currentFile))
        for di in dirs:
            shutil.rmtree(os.path.join(root,di))
        
        

def compare_protein_sequence(reference_seq,target_seq):
    """
    Function
    from check_variant_transcription if compared nt sequence are equal: Translate nt Bio.seq sequence, perform clustalw alignment and report aa changes between two sequences
    Input:
    - reference_seq['Bio.Seq.Seq'>]: reference nuclotide sequence
    - target_fasta_sequence['Bio.Seq.Seq'> ]: target nuclotide sequence
    Returns
    final_report[dict]: key[integer]= position value[str]= amino acid changes
    e.g.  {145: 'no_stop', 37: 'I-X'}
    """
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

   
def clustalw_align(reference_seq,target_seq,output):
    """
    Function
    from check_variant_transcription if compared nt sequence not equal: Translate nt Bio.seq sequence, perform clustalw alignment and report aa changes between two sequences
    Input:
    - reference_seq['Bio.Seq.Seq'>]: reference nuclotide sequence
    - target_fasta_sequence['Bio.Seq.Seq'> ]: target nuclotide sequence
    Returns
    final_report[dict]: key[integer]= position value[str]= amino acid changes
    e.g.  {145: 'no_stop', 37: 'I-X'}
    """

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

def check_variant_transcription(reference_seq,target_fasta_sequence,outdir):
    """
    Function
    from best_hit function: take nt sequence of reference and target seq, translate into proteins and compare aa sequence by checking elements in sequence lists or by clustalw alignment if different length 
    Input:  
    - reference_seq['Bio.Seq.Seq'>]: reference nuclotide sequence
    - target_fasta_sequence['Bio.Seq.Seq'> ]: target nuclotide sequence
    - outdir[str]: full path to the output dir  
    
    Returns
    mutations[dict]Position and change of mutated amino acid sequence including non stop or early stop
    e.g.  {145: 'no_stop', 37: 'I-X'}
    """
  
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

def pos_of_int_to_list(pos_of_int):
    """
    Function
    Take position of mutations of interest from workflow.txt file and return a list of position
    Input[str] = eg. 5,7,10-15,20
    return[list] = eg. [5,7,10,11,12,13,14,15,20]
    """
    if pos_of_int != 'N':
        total_list = []
        pos_int_list = pos_of_int.replace(" ","").strip("\n").split(",")
        for pos_dt in pos_int_list:
            region = pos_dt.split("-")
            if len(region) == 2 and int(region[0]) < int(region[1]):
                total_list.extend(range(int(region[0]),int(region[1])))
            elif len(region) == 1:
                total_list.append(int(pos_dt))
            else:
                print "format for mutations in workflow file is not correct eg. 5-20,46,75-90"
                sys.exit(1)
        return total_list
    else:
        return None

def best_hit(workflow_file,pileup_hash,path_to_reference_folder,outdir):
    """
    Function
    Take pileup_hash [dict] generated in extract_quality_metrics function and report only best hits if multiple reference sequences exist from a gene
    update the dictionary with info from workflow.txt eg. desrciption, report_type
    calculate control-coverage and mix_indicator values from chromosoml targets
    Input:
    - pileup_hash [dict] : generated in extract_quality_metrics function
    - path_to_reference_folder[str]: full path to the refrence folder 
    - outdir[str]: full path to the output dir 
    - workflow_file[str]: full path tp workflow.txt
    Returns
    -xml_output[dict]: pileup_hash dictionary updated with best_hits info for variants
    -control_coverage_value[float]: average depth from all chromosomal reference genes or 'failed' in chromosomal genes are not detected
    -mix_indicator[str] = nb of mix positions in chromosomal reference genes eg. (5:5052) (5 mix positions: total nb of chromosomal positions) or 'failed' in chromosomal genes are not detected
    """ 
    total_mix ={}
    control_coverage = []
    list_of_keys = pileup_hash.keys()
    xml_output = {}
    with open(workflow_file) as workflow_info:
        for line in workflow_info:
            final_output = []
            gene,ref_homology,workflow,pos_of_int,description,reporting_type,reporting_order = line.strip("\n").split("\t")            
            report_type = reporting_type + "_" + reporting_order
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
                    if total_mix.has_key('total_chromosomal_length'):
                        total_mix['total_chromosomal_length'].append(int(pileup_hash[best_hits[0]]['allele_length']))
                    else:
                        total_mix['total_chromosomal_length'] = [int(pileup_hash[best_hits[0]]['allele_length'])]
                    if total_mix.has_key('total_mix_positions'):
                        total_mix['total_mix_positions'].append(int(len(pileup_hash[best_hits[0]]['positions_mix'].keys())))
                    else:
                        total_mix['total_mix_positions'] = [int(len(pileup_hash[best_hits[0]]['positions_mix'].keys()))]
                    
                    mutation_info = {}
                    reference_seq = hash_reference_seq[best_hits[0]]
                    if pileup_hash[best_hits[0]]['coverage'] == float(100):
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
                            position_of_interest = pos_of_int_to_list(pos_of_int)
                            if position_of_interest != None:
                                relevent_position = dict([(key,mutation_info['positions_prot_modifications'][key]) for key in mutation_info['positions_prot_modifications'] if int(key) in position_of_interest])
                                pileup_hash[best_hits[0]]['relevent_position'] = relevent_position
                            else:
                                pileup_hash[best_hits[0]]['relevent_position'] = "ND"
                        else:
                            pileup_hash[best_hits[0]]['relevent_position'] = {}
                    else:
                        pileup_hash[best_hits[0]]['relevent_position'] = "ND"
                        
                    pileup_hash[best_hits[0]]['ref_cut_off'] = int(ref_homology)
                    pileup_hash[best_hits[0]]['report_type'] = str(report_type)
                    pileup_hash[best_hits[0]]['description'] = str(description)
                    xml_output[best_hits[0]] = pileup_hash[best_hits[0]]
                elif workflow == 'variant':
                    reference_seq = hash_reference_seq[best_hits[0]]
                    if pileup_hash[best_hits[0]]['coverage'] == float(100): #and best_hits[0].split("_")[0] != "23s":
                        alterations = check_variant_transcription(reference_seq,target_fasta_sequence,outdir)
                        pileup_hash[best_hits[0]]['protein_alterations'] = alterations
                    else:
                        pileup_hash[best_hits[0]]['protein_alterations'] = 'ND'
                    pileup_hash[best_hits[0]]['ref_cut_off'] = int(ref_homology)
                    pileup_hash[best_hits[0]]['report_type'] = str(report_type)
                    pileup_hash[best_hits[0]]['description'] = str(description)
                    xml_output[best_hits[0]] = pileup_hash[best_hits[0]]
                    
                    
                elif workflow == 'regulator':
                    pileup_hash[best_hits[0]]['nucleotide_alterations'] = {}
                    total_position_mismatch = pileup_hash[best_hits[0]]['position_nuc_mismatchs']
                    total_coverage_info = pileup_hash[best_hits[0]]['position_with_coverage']                    
                    position_of_interest = pos_of_int_to_list(pos_of_int)
                    if position_of_interest != None:
                        if position_of_interest != 0:
                            for pos_reg in position_of_interest:
                                if total_coverage_info.has_key(pos_reg):
                                    if total_position_mismatch.has_key(pos_reg):
                                        pileup_hash[best_hits[0]]['nucleotide_alterations'][pos_reg] = total_position_mismatch[pos_reg]
                                else:
                                    pileup_hash[best_hits[0]]['nucleotide_alterations'][pos_reg] = 'ND'
                    else:
                        pileup_hash[best_hits[0]]['nucleotide_alterations'] = 'ND'
                    pileup_hash[best_hits[0]]['ref_cut_off'] = int(ref_homology)
                    pileup_hash[best_hits[0]]['report_type'] = str(report_type)
                    pileup_hash[best_hits[0]]['description'] = str(description)
                    xml_output[best_hits[0]] = pileup_hash[best_hits[0]]
            else:
                if workflow == 'mutant':
                    xml_output[gene] = {}
                    xml_output[gene]['housekeeping_detection'] = 'ND'
                    report_type = reporting_type + "_" + reporting_order
                    xml_output[gene]['report_type'] = str(report_type)
                    xml_output[gene]['description'] = str(description)              
                elif reporting_type == 'required' and workflow != 'mutant':
                    xml_output[gene] = {}
                    xml_output[gene]['absolute_reporting'] = 'ND'
                    report_type = reporting_type + "_" + reporting_order
                    xml_output[gene]['report_type'] = str(report_type)
                    xml_output[gene]['description'] = str(description)
                else:
                    pass
            
    if len(control_coverage) > 0:       
        control_coverage_value = sorted(control_coverage)[0]
    else:
        control_coverage_value = 'failed'
    if len(total_mix) > 0:
        mix_indicator = str(sum(total_mix['total_mix_positions'])) + ":" + str(sum(total_mix['total_chromosomal_length']))
    else:
        mix_indicator = 'failed'    
    for key in xml_output.keys():
        for sub_key in xml_output[key].keys():
            if sub_key == 'position_insertions' or sub_key == 'position_deletions':
                print sub_key
                print xml_output[key][sub_key]
            
    return xml_output,control_coverage_value,mix_indicator


def report_xml(outdir,workflow_name,version,prefix,xml_output,control_coverage_value,mix_indicator):
    """
    Function
    generate xml file from dictionary/info generated by best_hit function
    Input:  
    -outdir[str]: full path to the output dir
    -workflow_name[str]: species_workflow (eg. staphylococcus_typing
    -version[str]: version number
    -prefix[str]: sampleid (e.g. NGSLIMS_molis number, 9188_H14276008707)
    -xml_output[dict]:primary key=gene_id, value= dictionary with keys= type of metric and value= value metric
    -control_coverage_value[float]: from best_hit function
    
    - Returns
     generate and print xml file in output folder
    
    """ 

    genome_id = prefix.split(".")[0]   
    xml_log_file = open(outdir + "/" + genome_id + ".results.xml", "w")
    root = etree.Element("ngs_sample", id = genome_id)
    workflow = etree.SubElement(root, "workflow", value=workflow_name, version = version)
    results = etree.SubElement(root, 'results')
    #gene_profiling = etree.SubElement(results, 'gene_profiling')
    coverage_control = etree.SubElement(results, "coverage_control", value=str(control_coverage_value))
    mix_indicator = etree.SubElement(results, "mix_indicator", value=mix_indicator)   
    for gene_cluster, value in xml_output.items():
        #gene = gene_cluster.split('_')[0]
        gene = gene_cluster
        # print xml_output['position_deletions']
        # print xml_output['position_insertions']
        result = etree.SubElement(results, "result", type="gene", value = gene)
        if value.has_key('housekeeping_detection'):
            etree.SubElement(result, "result_data", type="mode", value = 'mutant')
            etree.SubElement(result, "result_data", type="detection", value = 'ND') #change here to replace output to not_detected
            etree.SubElement(result, "result_data", type="description", value = str(value['description']))
            etree.SubElement(result, "result_data", type="report_type", value = str(value['report_type']))
        elif value.has_key('absolute_reporting'):
            etree.SubElement(result, "result_data", type="mode", value = 'variant')
            etree.SubElement(result, "result_data", type="detection", value = 'ND') #change here to replace output to not_detected
            etree.SubElement(result, "result_data", type="description", value = str(value['description']))
            etree.SubElement(result, "result_data", type="report_type", value = str(value['report_type']))        
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
                etree.SubElement(result, "result_data", type="detection", value = 'D')# change here to replace output for detection
                gene_predicted_seq = []
                for contig_predicted in value['predicted_seq']:
                    gene_predicted_seq.append(str(contig_predicted.id))
                    gene_predicted_seq.append(str(contig_predicted.seq))
                fasta_contigs_file = open(outdir + '/' + gene + '_contigs.fasta',"w")
                fasta_contigs_file.write('\n'.join(gene_predicted_seq))
                fasta_contigs_file.close
                                
            elif value['homology'] >= float((value['ref_cut_off'])) and value['coverage'] != float(100):# can relaxe the condition of uncertain here
                etree.SubElement(result, "result_data", type="detection", value = 'U') # change here to moditify output for uncertain
                all_contigs = []
                for contig_predicted in value['predicted_seq']:
                    all_contigs.append(str(contig_predicted.id))
                    all_contigs.append(str(contig_predicted.seq))
                fasta_contigs_file = open(outdir + '/' + gene + '_contigs.fasta',"w")
                fasta_contigs_file.write('\n'.join(all_contigs))
                fasta_contigs_file.close 
            else:
                etree.SubElement(result, "result_data", type="detection", value = 'ND')# change here to modidy output for absent
            etree.SubElement(result, "result_data", type="description", value = str(value['description']))
            etree.SubElement(result, "result_data", type="report_type", value = str(value['report_type']))
            etree.SubElement(result, "result_data", type="coverage", value = str(value['coverage']))
            etree.SubElement(result, "result_data", type="homology", value = str(value['homology']))
            etree.SubElement(result, "result_data", type="depth", value = str(value['depth']))
            etree.SubElement(result, "result_data", type="coverage_distribution", value = str(value['coverage_distribution']))
            #XXXXXXX
            if value.has_key('position_insertions'):
                insertion_pos =[]
                if len(value['position_insertions']):
                    for position,mutation in sorted(value['position_insertions'].items()):
                        relevent_info =  str(position) + ':' + mutation
                        insertion_pos.append(relevent_info)
                insertion_pos_report = ';'.join(insertion_pos)
                if len(insertion_pos_report) > 0:
                    etree.SubElement(result, "result_data", type="insertions", value = insertion_pos_report)
                else:
                    etree.SubElement(result, "result_data", type="insertions", value = "None")
            if value.has_key('position_deletions'):
                deletion_pos =[]
                if len(value['position_deletions']):
                    for position,mutation in sorted(value['position_deletions'].items()):
                        relevent_info =  str(position) + ':' + mutation
                        deletion_pos.append(relevent_info)
                deletion_pos_report = ';'.join(deletion_pos)
                if len(deletion_pos_report):
                    etree.SubElement(result, "result_data", type="deletions", value = deletion_pos_report)
                else:
                    etree.SubElement(result, "result_data", type="deletions", value = "None")
            #XXXXXXX
            if len(value['positions_mix']) > 0:
                report = []
                for position,mutation in sorted(value['positions_mix'].items()):
                    all_values = mutation.values()
                    new = [str(elem) for elem in all_values]
                    mismatch_info = str(position) + ':' + '-'.join(mutation.keys()) + '-' + '-'.join(new)
                    report.append(mismatch_info)
                mix = ';'.join(report)
                etree.SubElement(result, "result_data", type="mix", value = mix)
            else:
                etree.SubElement(result, "result_data", type="mix", value = 'None')
            
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
    print etree.tostring(root, pretty_print=True)
    print >> xml_log_file, etree.tostring(root, pretty_print=True)
    
