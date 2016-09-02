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
	       - bowtie2/2.1.0patch
	       - numpy/python2.7/1.7.1
	       - pysam/python2.7/0.7.5
	       - phe/gene_finder/1-0
	       - samtools/0.1.18
	       - lxml/python2.7/3.2.3
	       - yaml/1.1
	       

    usage[1] has two required arguments
      1) input_directory (-i) containing a pair of processed fastq files (format *.processed.*.fastq).
      2) workflow (-w) e.g 'staphylococcus-aureus-typing.ngsservice'
    usage[1] optional arguments include
      3) full path to reference dir containing the reference genbank and fasta files as reference.fasta(multifasta), 
        workflow.txt(tab_delimited). This is set by default to the latest reference dir for the active workflow when 
        module phe/gene_finder/<major-minor>prereq is loaded; otherwise a path needs to be passed to the 
        -gf|--gene_file_directory param.
	    4) output_directory (-o) is optional; the default creates an output dir within the input directory

	
	usage[2] has two required arguments
	  1) input_directory should contain only two paired fastq files (format *.fastq or *.fastq.gz)
      2) full path to reference dir containing the reference genbank and fasta files as reference.fasta(multifasta), 
        workflow.txt(tab_delimited) needs to be passed to the -gf|--gene_file_directory param.


	usage[3] has three required arguments
	  1) full path to fastq_1, preferably processed, passed to -1|-fastq_1 param
	  2) full path to fastq_2, preferably processed, passed to -2|-fastq_2 param
      3) full path to reference dir containing the reference genbank and fasta files as reference.fasta(multifasta), 
        workflow.txt(tab_delimited) needs to be passed to the -gf|--gene_file_directory param.
	
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
        print file_description + " detected;" + filepath

parser = argparse.ArgumentParser()
parser.add_argument('-input_directory', '-i', help='please provide an input directory')
parser.add_argument('-workflow', '-w', help='If using a workflow you must specify an input directory')
parser.add_argument('-fastq_1', '-1', help='Fastq file pair 1')
parser.add_argument('-fastq_2', '-2', help='Fastq file pair 2')
parser.add_argument('-output_directory', '-o',type=str,help='Output root directory')
parser.add_argument('-obn',
                    '-output_dir_basename',
                    help='pass serotyping component name; used to name output dir; over-ridden by any -output_dir|-o arg; \n[Intended for use in automated qsub_script calls to workflow entry point, but implemented for all three entry points]')
parser.add_argument('--gene_file_directory',
                    '-gf',
                    default=os.environ.get('GENE_FINDER_REFERENCE_DIR'),
                    #tmp#default=os.environ.get('ESCHERICHIA_COLI_SEROTYPING_REFERENCE_DIR'),
                    help='Path to dir containing a multientry reference.fasta file and a workflow.txt file.\nSet automatically by gene_finder modulefile if env_var WORKFLOW_NAME has been set.')
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

        # GENE_FINDER_REFERENCE_DIR is an environment variable assigned by module 
        # /phe/phe_ref/<workflow_name>/gene_finder_reference
        # This phe_ref module is loaded by the phe/gene_finder module, given the 
        # pre-requisite assignment of the workflow_name (i.e. workflow.split(.)[0]) 
        # to the environment variable WORKFLOW_NAME as below
        # export WORKFLOW_NAME=workflow_name

        ## If the pre-requisite is not met, the environment variable 
        ## GENE_FINDER_REFERENCE_DIR is not set, leaving opt.gene_file_directory == None. 
        ## That condition needs amended before proceeding
        if not opts.gene_file_directory:
            if not opts.obn:
                print("If you are using the 'workflow' entry point but have not loaded gene_finder using the 'prereq' module, you need also to pass to the -gf|--gene_file_directory param a path to a reference_dir containing files reference.fasta and workflow.txt")
                sys.exit(1)
            else:
                # output_dir_basename is component name; use to map to ref_dir
                opts.gene_file_directory = os.environ[
                                       '{}_REFERENCE_DIR'.format(opts.obn.upper())] # tracetrace if WORKFLOW_NAME has not been exported before modulefile is loaded

        if len(fastq_files) < 2:
            print "Fastq files are not pairs!"
            sys.exit(1)
        if not opts.output_directory:
            if not opts.obn:
                opts.output_directory = opts.input_directory + '/gene_finder'
            else:
                opts.output_directory = '{}/{}'.format(opts.input_directory, opts.obn)
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
        print 'opts:',opts #tmp#


    # second option: need to provide forward and reverse fastq_files using -1 and -2 options
    # and path to reference dir contining reference.fasta and workflow.txt files
    elif opts.fastq_1 or opts.fastq_2:
        check_file_exists(opts.fastq_1, 'Fastq 1')
        check_file_exists(opts.fastq_2, 'Fastq 2')
        if not opts.gene_file_directory:
            if not opts.obn:
                print("If you are passing full processed fastq paths to the fastq-1 and fastq-2 params, you need also to pass to the -gf|--gene_file_directory param a path to a reference_dir containing files reference.fasta and workflow.txt")
                sys.exit(1)
            else:
                # output_dir_basename is component name; use to map to ref_dir
                opts.gene_file_directory = os.environ[
                                       '{}_REFERENCE_DIR'.format(opts.obn.upper())] # tracetrace if WORKFLOW_NAME has not been exported before modulefile is loaded
        check_file_exists(opts.gene_file_directory, 'gene_file_directory')

        if not opts.output_directory:
            if not opts.obn:
                opts.output_directory = os.path.dirname(opts.fastq_1) + '/gene_finder'
            else:
                opts.output_directory = '{}/{}'.format(os.path.dirname(opts.fastq_1), opts.obn)
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
        print 'opts:',opts #tmp#

    # third option:need to provide the path for a dir containing the paired fastq files 
    elif opts.input_directory:
        check_file_exists(opts.input_directory, 'input_directory')
        if not opts.gene_file_directory:
            if not opts.obn:
                print("If you are passing a directory containing processed fastqs as arg to the -i|--input_directory param, you need also to pass to the -gf|--gene_file_directory param a path to a reference_dir containing files reference.fasta and workflow.txt")
                sys.exit(1)
            else:
                # output_dir_basename is component name; use to map to ref_dir
                opts.gene_file_directory = os.environ[
                                       '{}_REFERENCE_DIR'.format(opts.obn.upper())] # tracetrace if WORKFLOW_NAME has not been exported before modulefile is loaded
        check_file_exists(opts.gene_file_directory, 'gene_file_directory')

        # glob preferentially for processed fastq files
        fastq_files = glob.glob('%s/%s' % (opts.input_directory,
                                           glob_pattern))
        # if no processed fastq then use whatever's available
        if not fastq_files:
            fastq_files = glob.glob('%s/%s' % (opts.input_directory,
                                               "*fastq*"))
        if not opts.output_directory:
            if not opts.obn:
                opts.output_directory = opts.input_directory + '/gene_finder'
            else:
                opts.output_directory = '{}/{}'.format(opts.input_directory, opts.obn)
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
        print 'opts:',opts #tmp#


    fasta_file = opts.gene_file_directory + '/reference.fasta'
    workflow_file= opts.gene_file_directory + '/workflow.txt'
    for file in [ fasta_file, workflow_file ]:
       check_file_exists(file, os.path.basename(file) )


    outdir = opts.output_directory
    cut_off = opts.cut_off
    bowtie_options = opts.bowtie_options
    minimum_coverage = opts.min_cov
    component = 'gene_finder' if not opts.obn else opts.obn
    stderr_log_output = opts.log_directory + "/" + component + ".stderr"
    stdout_log_output = opts.log_directory + "/" + component + ".stdout"
    logger = log_writer.setup_logger(stdout_log_output, stderr_log_output)

    if workflow_name != None:
        gene_finder_functions.run_gene_finder(logger,
                                    outdir,
                                    fasta_file,
                                    workflow_file,
                                    fastq_files,
                                    bowtie_options,
                                    ids,
                                    cut_off,
                                    minimum_coverage,
                                    opts.log_directory,
                                    workflow_name=workflow_name,
                                    version=version)
#        try_and_except(stderr_log_output,
#                       write_component_complete,
#                       opts.output_directory)
        return 0

    elif workflow_name == None:
        print "this condition"
        gene_finder_functions.run_gene_finder(logger,
                                    outdir,
                                    fasta_file,
                                    workflow_file,
                                    fastq_files,
                                    bowtie_options,
                                    ids,
                                    cut_off,
                                    minimum_coverage,
                                    opts.log_directory,
                                    workflow_name="None",
                                    version="None")

if __name__ == '__main__':
    rc = main()
    if rc is None:
        pass
    elif rc:
        exit(rc)
    else:
        write_component_complete(opts.output_directory)

