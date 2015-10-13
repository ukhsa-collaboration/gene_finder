#! /usr/bin/env python

def help_function():
	print """
	**gene finder**
	================================================================================
	command: python gene_finder.py [usage]
	    [usage1]: -w [workflow] -i [input_directory] -o [output_directory:default in input directory]
	    [usage2]: -i [input_directory] -gf [path to gene reference directory] -o [output_directory:default in input directory]
	    [usage3]: -1 [path fastq 1] -2 [path fastq 2] -gf [path to gene reference directory] -o [output_directory:default in fastq directory]
	
	       optional
	       -bowtie-options:(list),default=['-q','--very-sensitive-local','--no-unal','-a']
	       -cut_off':(str), cut-off used to identify mix (first integer) and indels (second integer), default='86:50' mix is considered if non-identical reads are < 86% of total reads, indel is considered if  > 50% of total reads indicate so.
	       -min_cov:(int), cut_off used to define minimum acceptable coverage at each position, default=5 acceptable coverage correpond to 5% of maximum coverage detected among all positions for a given gene.
	       
	       module dependancies:
	       - python/2.7.5
	       - clustalw/2.1
	       - biopython/python2.7/1.61
	       - phe/common_modules/1-0
	       - bowtie2/2.1.0
	       - numpy/python2.7/1.7.1
	       - pysam/python2.7/0.7.5
	       - phe/gene_finder/1-0
	       - samtools/0.1.18
	       - lxml/python2.7/3.2.3
	       - yaml/1.1
	       
    usage[1] need at least two arguments
        1) input_directory (-i) containing a pair of processed fastq files (format *.processed.*.fastq).
	    2) workflow (-w) e.g 'staphylococcus-aureus-typing-1.0' is folder [reference.fasta(multifasta), workflow.txt(tab_delimited)
	         the config file for the work flow is hard coded but can be overridden using the -c argument
	    3) output_directory (-o) is optional, default an output file is created in the input directory
	
	usage[2]
	    1) input_directory should contain only two paired fastq files (format *.fastq or *.fastq.gz)
	
	Reference directory should contain two files:
	 - multifasta file named "reference.fasta" containing nucleotide sequences with gene name in header: eg. gyrA or sea_1, sea_2 for report best match
	 - a  tab delimited file named "workflow.txt" containing
	     gene_id : string (eg. gyrA,aac(6')-Ib-cr), or sea if multiple sequneces for the same gene are in reference.fasta(eg. sea_1, sea_2)
	     percent-homology: cut off used to determine presence absence (eg. 90)
	     mode: three mode are accepted variant, mutant, regulator
	     mutation-position : positions be reported as integers seperated by (,), consecutive positions are separated by (-) (eg. 50,60-106,70,75...)
	     description:    string gene description (eg. resistance, toxin ...)
	     report_type: string "required" to report presence, absence or uncertain or "conditional" to report only presence or uncertain
	     order of reporting: intger 1,2 ....
	
	"""
	
# importing system libraries
import os, os.path, sys, subprocess, getopt, argparse, glob, yaml, inspect

module_folder_paths = ["modules"]
for module_folder_path in module_folder_paths:
	module_folder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],module_folder_path)))
	if module_folder not in sys.path:
		sys.path.insert(1, module_folder)
from utility_functions import *
import log_writer
import gene_finder_functions
import generate_mpileup_file
import parsing_mpileup
import extract_quality_metrics

def check_file_exists(filepath, file_description):
        if not os.path.exists(filepath):
            print(file_description + " (" + filepath + ") does not exist")
            sys.exit(1)
        else:
            print file_description + " detected"
		

parser = argparse.ArgumentParser()
parser.add_argument('-input_directory', '-i', help='please provide an input directory')
parser.add_argument('-workflow', '-w', help='If using a workflow you must specify an input directory')
parser.add_argument('-fastq_1', '-1', help='Fastq file pair 1')
parser.add_argument('-fastq_2', '-2', help='Fastq file pair 2')
parser.add_argument('-output_directory', '-o',type=str,help='Output root directory')
#parser.add_argument('-gene_file_directory', '-gf',help='Directory containing the reference genbank and fasta files')
parser.add_argument('-gene_file_directory',
                    '-gf',
                    default=os.environ.get('GENE_FINDER_REFERENCE_DIR'),
                    help='Directory containing the reference genbank and fasta files')
parser.add_argument('-bowtie_options','-bopt',help='path_to_clustalw',default=['-q','--very-sensitive-local','--no-unal','-a'])
parser.add_argument('-cut_off', '-c',type=str, help='cut_off values Sub:Ind',default='84:50')
parser.add_argument('-min_cov', '-m',type=str, help='min nb of reads for coverage confidence',default='5')
parser.add_argument('--log_directory', '-log', help='please provide the path for log directory')

opts = parser.parse_args()

def main():
    if (len(sys.argv) <= 1):
        parser.print_help()
        sys.exit()
    else:
	fastq_files = []	
	glob_pattern = "*.processed*.fastq*"
	ids = None
	workflow_name = None
	version = None

    if opts.workflow:
        if not opts.input_directory:
            print("If using a workflow you must specify an input directory")
            sys.exit(1)
        else:
            check_file_exists(opts.input_directory, 'input directory')
            fastq_files = glob.glob(opts.input_directory + "/" + glob_pattern)
        #opts.gene_file_directory = os.path.dirname(os.path.realpath(__file__)) + "/gene_reference_data/" + opts.workflow

        # GENE_FINDER_REFERENCE_DIR is assigned by module /phe/phe_ref/<workflow_name>/gene_finder_reference/module_version
        # loading this phe_ref module by the phe/gene_finder module requires the prior command-line assignment of the 
        # workflow_name (i.e. workflow.split(.)[0]) to the environment variable WORKFLOW_NAME as below
        # export WORKFLOW_NAME=workflow_name

        ## 2015OCT13 now set as default within parser
        ## opts item will be None if os.environ is not set
        ## opts.gene_file_directory = os.environ['GENE_FINDER_REFERENCE_DIR'] # see preceeding lines if KeyError: 'GENE_FINDER_REFERENCE_DIR'

        if len(fastq_files) < 2:
            print "Fastq files are not pairs!"
            sys.exit(1)
        if not opts.output_directory:
            opts.output_directory = opts.input_directory + '/gene_finder'
            if not os.path.isdir(opts.output_directory): os.makedirs(opts.output_directory) #make output_directory

        (seqDir,seqFileName) = os.path.split(fastq_files[0])
        split_seqFileName = seqFileName.split('.')
        version = split_seqFileName[-5]
        workflow_name = split_seqFileName[-6]
        ids = seqFileName.split(workflow_name)[0]
        ids = ids[:-1]
        if not opts.log_directory:
            opts.log_directory = opts.input_directory + '/logs'
            if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directoryectory


    # second option: need to provide forward and reverse fastq_files using -1 and -2 options
    elif opts.fastq_1 or opts.fastq_2:
        check_file_exists(opts.fastq_1, 'Fastq 1')
        check_file_exists(opts.fastq_2, 'Fastq 2')
        if not opts.gene_file_directory:
            print("If you are using -1 and -2 options, you must provide the reference fasta file and GF workflow directory folder")
            sys.exit(1)
        else:
            check_file_exists(opts.gene_file_directory, 'gene_file_directory')

        if not opts.output_directory:
            opts.output_directory = os.path.dirname(opts.fastq_1) + '/gene_finder'
            if not os.path.isdir(opts.output_directory):
                os.makedirs(opts.output_directory)
                opts.log_directory = os.path.dirname(opts.fastq_1) + '/logs'
                if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directory  
                
        else:
            if not os.path.isdir(opts.output_directory):
                os.makedirs(opts.output_directory)
                opts.log_directory = opts.output_directory + '/logs'
                if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directory
                
        fastq_files.append(opts.fastq_1)
        fastq_files.append(opts.fastq_2)
        gene_file_directory = opts.gene_file_directory
        (SeqDir,seqFileName) = os.path.split(fastq_files[0])
        (ids,ext) = os.path.splitext(seqFileName)

    # third option:need to provide the path for a dir containing the paired fastq files 
    elif opts.input_directory:
        check_file_exists(opts.input_directory, 'input_directory')
        if not opts.gene_file_directory:
            print("If you are providing the input directory for the fastqs, the reference fasta file and gene reference folder which should contain reference.fasta and workflow.txt files")
            sys.exit(1)
        else:
            check_file_exists(opts.gene_file_directory, 'gene_file_directory')
        fastq_files = glob.glob(opts.input_directory + "/*fastq*")
        if not opts.output_directory:
            opts.output_directory = opts.input_directory + '/gene_finder'
            if not os.path.isdir(opts.output_directory): os.makedirs(opts.output_directory) #make output_directory
            opts.log_directory = opts.input_directory + '/logs'
            if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directory  
        else:
            os.makedirs(opts.output_directory)
            opts.log_directory = opts.output_directory + '/logs'
            if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directory
        
        (seqDir,seqFileName) = os.path.split(fastq_files[0])
        (ids,ext) = os.path.splitext(seqFileName)
        if not opts.log_directory:
            opts.log_directory = opts.output_directory + '/logs'
        if not os.path.isdir(opts.log_directory): os.makedirs(opts.log_directory) #make log_directory
	    

    fasta_file = opts.gene_file_directory + '/reference.fasta'
    workflow_file= opts.gene_file_directory + '/workflow.txt'
    for file in [ fasta_file, workflow_file ]:
       check_file_exists(file, os.path.basename(file) )
 
        
    outdir = opts.output_directory
    cut_off = opts.cut_off
    bowtie_options = opts.bowtie_options
    minimum_coverage = opts.min_cov
    stderr_log_output = opts.log_directory + "/" + 'gene_finder'+ ".stderr"
    stdout_log_output = opts.log_directory + "/" + 'gene_finder'+ ".stdout"
    logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)
    if workflow_name != None:
        gene_finder_functions.run_gene_finder(logger,outdir,fasta_file,workflow_file,fastq_files,bowtie_options,ids,cut_off,minimum_coverage,opts.log_directory,workflow_name=workflow_name,version=version)
        try_and_except(stderr_log_output,write_component_complete,opts.output_directory)
    elif workflow_name == None:
        print "this condition"
        gene_finder_functions.run_gene_finder(logger,outdir,fasta_file,workflow_file,fastq_files,bowtie_options,ids,cut_off,minimum_coverage,opts.log_directory,workflow_name="None",version="None")

main()

