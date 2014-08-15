#!/usr/bin/env python
import logging
from subprocess import call
import os,sys,inspect,re,subprocess

module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
	module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
	if module_folder not in sys.path:
		sys.path.insert(1, module_folder)
import utility_functions
import log_writer

def index_reference(fasta_file,logger):
    #==================================================================================#
    ## Description: Check and Build Bowtie2_index_files if not yet created##
    ## Input: 1- string: full path to reference fasta file
    ##        2- logger: defined as "log_writer.setup_logger(stdout_log_output, stderr_log_output)"
    ## Output: Index files in the fasta_file directory
    #==================================================================================#
    bowtie2_index = fasta_file + '.1.bt2'
    if os.path.exists(bowtie2_index):
        log_writer.info_header(logger,'bowtie2 reference index files already exist')
    else:
        process = subprocess.Popen(['bowtie2-build', '-f', fasta_file, fasta_file], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        process.wait()
        log_writer.log_process(logger, process, log_error_to = "bowtie2: reference indexing successfully completed")
        log_writer.info_header(logger,'bowtie2: reference indexing successfully completed')
        
def samtools_faidx(fasta_file,logger):
    #==================================================================================#  
    ## Description: Check and Create SAMtools index files if not yet created##
    ## Input: 1- string => full path to reference fasta file
    ##        2- logger: defined as "log_writer.setup_logger(stdout_log_output, stderr_log_output)"
    ## Output: Index files (.fai) in the fasta_file directory
    #==================================================================================#
    fai_index = fasta_file + '.fai' 
    if os.path.exists(fai_index):
        if os.path.getsize(fai_index) > 0:
            log_writer.info_header(logger,'samtools reference fasta faidx file exist')
        else:
            print('samtootls faidx output is empty, check reference fasta format')
            log_writer.info_header(logger,'samtootls faidx output is empty, check reference fasta format')
            sys.exit(1)
    else:
        process = subprocess.Popen(['samtools', 'faidx', fasta_file], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        process.wait()
        log_writer.log_process(logger, process, log_error_to = "info")
        if os.path.getsize(fai_index) > 0:
            print("samtootls faidx sucsessfully created")
            log_writer.info_header(logger,'samtootls faidx sucsessfully created')
        else:
            print('samtootls faidx output is empty, check reference fasta format')
            log_writer.info_header(logger,'samtootls faidx output is empty, check reference fasta format')
            sys.exit(1) 

def run_bowtie_on_indices(fasta_file,forward_fastq,reverse_fastq,outdir,workflow_name,version,prefix,bowtie_options,logger):
    #====================================================================================================# 
    ## Description: Create tmp file to write outputs, Run bowtie2 with default_options =['-q','--very-sensitive-local','--no-unal','-a']
    ## Inputs: -1 string => path_to 1)reference_fasta 2) paired_fastq_forward 3) paired_fastq_reverse 4) output_directory
    ##         -2 string => species_workflow (eg. staphylococcus_typing), version (eg. 1-0), prefix (eg. molis id) used for labelling all_outputs
    ##         -3 list bowtie_options (if not default)
    ##         -4 logger defined as "log_writer.setup_logger(stdout_log_output, stderr_log_output)"
    ## Output -1 SAM file in tmp file
    ## NB: bowtie_options can be modified using argument --bopt (-bowtie_options) with options provided as a list##
    #====================================================================================================#  
    sam_output = outdir + "/tmp/" + prefix + '.sam'
    #if os.path.exists(sam_output):
    #    print "SAM file already created for this sample ??"
    #    log_writer.info_header(logger,"SAM file already created for this sample ??")
    #    sys.exit(1) 
    #else:
    command = ['bowtie2']
    command += ['-1', forward_fastq, '-2', reverse_fastq,'-x',fasta_file,'-S', sam_output]
    command += bowtie_options
    log_writer.info_header(logger,'run bowtie command')
    log_writer.info_header(logger, ' '.join(map(str,command)))
    process = subprocess.Popen(command,stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    process.wait()
    log_writer.info_header(logger,'SAM file successfully created')

def modify_bowtie_sam(path_sam_folder,prefix,logger):
    #==================================================================================# 
    ## Description: Deduce 256 from the second column to unset secondary alignment. Reads with 250 bit scrore are not reported when using Samtools pileup
    ## Input: -1 string => path_sam_folder
    ##        -2 string => prefix
    ##        -3 logger defined as "log_writer.setup_logger(stdout_log_output, stderr_log_output)"
    ## Return sam.mod file in the same directoruy of the original .sam file
    #==================================================================================#    
    path_sam_file = path_sam_folder + "/" + prefix + ".sam"
    with open(path_sam_file) as sam, open(path_sam_file + '.mod', 'w') as sam_mod:
        for line in sam:
            if not line.startswith('@'):
                fields = line.split('\t')
                flag = int(fields[1])
                flag = (flag - 256) if (flag > 256) else flag
                sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
            else:
                sam_mod.write(line)
    log_writer.info_header(logger,'SAM file successfully mofified to unset secondary alignment')

def run_samtools_bam(path_to_tmp_file,fasta_file,prefix,logger):
    #==================================================================================# 
    ## Description: Run SAMtools and generate BAM, sort and index BAM files
    ## Input: -1 string => path_to_tmp_file containing SAM file, string : prefix for naming_SAMtools_outputs
    ##        -2 string prefix 
    ##        -3 logger defined as "log_writer.setup_logger(stdout_log_output, stderr_log_output)"
    ## output : generate BAM, sorted BAM and pileup in tmp_file
    #==================================================================================#  
    sam_output = path_to_tmp_file + "/"+ prefix + ".sam"
    bam_output = path_to_tmp_file + "/" + prefix + ".bam"
    bam_sorted_output = path_to_tmp_file + "/" + prefix + ".sorted"
    bam_sorted_index_output = bam_sorted_output + ".bam"
    pileup_output = path_to_tmp_file + "/" + prefix + '.pileup'
        
    #if os.path.exists(path_to_tmp_file + "/" + prefix + ".pileup"):
    #    print "pileup file already created for this sample ??"
    #    log_writer.info_header(logger,"SAM file already created for this sample ??")
    #    sys.exit(1)
    #else:
            
    log_writer.info_header(logger,'converting .sam to .bam')
    process = subprocess.Popen(['samtools', 'view', '-b', '-o', bam_output, '-q', '1', '-S', sam_output + '.mod'],stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    log_writer.info_header(logger,'samtools sorting .bam file')
    process = subprocess.Popen(['samtools', 'sort', bam_output, bam_sorted_output],stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    log_writer.info_header(logger,'samtools indexing .bam file')
    process = subprocess.Popen(['samtools', 'index', bam_sorted_index_output],stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    process.wait()
    log_writer.log_process(logger, process, log_error_to = "info")
    log_writer.info_header(logger,'generating mpileup file')
    log_writer.info_header(logger,' '.join(map(str,['samtools', 'mpileup', '-A', '-B', '-f', fasta_file, bam_sorted_index_output])))
    with open(pileup_output, 'w') as sam_pileup:
        process = subprocess.Popen(['samtools', 'mpileup', '-A', '-B', '-f', fasta_file, bam_sorted_index_output], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        for x in process.stdout:
            sam_pileup.write(x)
        process.wait()
        log_writer.log_process(logger, process, log_error_to = "couldn't generate mpileup file")


def generate_mpileup(fasta_file,forward_fastq,reverse_fastq,outdir,workflow_name,version,ids,bowtie_options,logger):
    #==================================================================================# 
    ## wrapper script with all steps needed to generate mpileup: generating tmp file for outputs, indexing reference, generating SAM, BAMs and mpileup
    #==================================================================================# 
    if os.path.exists(outdir + '/tmp'):
        log_writer.info_header(logger,"out_dir temp file already exist")
    else:
        log_writer.info_header(logger,"creating tmp file in the output directory")
        os.mkdir(outdir + '/tmp')
    stderr_log_output = outdir + "/" + ids+ ".stderr.log"
    stdout_log_output = outdir + "/" + ids+ ".stdout.log"
    logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)
    path_to_tmp_file = outdir + '/tmp'    
    run_bowtie_on_indices(fasta_file,forward_fastq,reverse_fastq,outdir,workflow_name,version,ids,bowtie_options,logger)
    modify_bowtie_sam(path_to_tmp_file,ids,logger)
    run_samtools_bam(path_to_tmp_file,fasta_file,ids,logger)
    
    
