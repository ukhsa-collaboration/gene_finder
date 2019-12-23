# GeneFinder

<h3>Dependencies</h3>

GeneFinder is written in Python 2.7.5 and requires the following dependencies:
- lxml 3.2.3, 
- numpy 1.7.1
- yaml 1.1
- pysam 0.7.5
- Biopython 1.61
- Samtools 1.0.18
- Bowtie 2.1.0
- Clustalw 2.1


<h3>Preparing reference sequences</h3>

The nucleotide sequences of genes of interest need to be stored in a multi FASTA format with only the names of genes in headers and MUST be designated ‘reference.fasta’. Gene names must NOT include underscores [‘_’]; these are reserved to designate variant sequences of a given gene (e.g gyrA_1, gyrA_2) which are treated by gene_finder as a sub-group with only the best matching variant selected for further processing and reporting. Before running gene_finder, an additional tab delimited table need to be prepared. This MUST be named ‘workflow.txt’ and should include for each gene sought seven description fields (see example below). Both files need to be located in the same folder which will be referred to in this document by ‘gene file directory’.


| Gene ID | Similarity Cutoff  | Mode  | Positions of Interest  | Description | Report Type | Reporting Posisiton
|---|---|---|---|---|---|---|
|  sea | 90 | variant  | N  | Enterotoxins  | required | 1
|  blaZ | 90 | variant  | N  | resistance  | conditional | 2
|  rpoB | 90 | mutant  | 455-566;575;700-760  | resistance  | required | 3
|  23s | 95 | regulator  | 2601-2605  | resistance  | required | 4


<h3>Running GeneFinder</h3>

gene_finder.py –1 <fastq_forward> -2 <fastq_reverse> -gf <reference_file_directory> -o <output_directory>

- -1	<string> full path to Illumina paired fastq forward
- -2	<string> full path to Illumina paired fastq reverse
- -o	<string> full path to output directory
- -gf	<string> full path to gene file directory which is expect to include two files called ‘reference.fasta’ and ‘workflow.txt’ prepared as described above

<h3>Outputs</h3>

GeneFinder output results for each sample in an XML format 

Gene_finder XML result includes two global control metrics:
**coverage_control**
Reports the average depth of reads for chromosomal genes, calculated from all genes labeled by mode ‘mutant’ in ‘workflow.txt’, return 'failed’ if no chromosomal gene is detected or included in reference sequences
Example value="122.68"
Reports the estimated average of read depth expected for chromosomal genes, thus present in one copy in the genome

**mix_indicator**
Reports the number of positions with mixture of alleles and total number of chromosomal positions inspected with values separated by colon, return 'failed' if chromosomal genes are not detected or included in the reference sequences
Example value="4:13928"
Example reporting 4 positions with mixture of nucleotides amongst 13928 positions representing the total number of positions of all chromosomal genes inspected.

For each gene sought, thereby specified in workflow.txt, gene_finder reports the following additional metrics:
**gene**
Returns gene name as defined in the FASTA reference header

**mode**
Returns mode of analysis (e.g variant, mutant or regulator) as defined in the workflow.txt

**detection**
Reports gene detection with three designations:
- ‘D’ for detected when coverage is 100 % and nucelotide similarity is greater than similarity cutoff value defined for the gene in the ‘workflow.txt’ file.
- ‘ND’ for not detected if similarity is below similarity cutoff value
- ‘U’ for uncertain if similarity is greater than similarity cutoff value but coverage is below 100%

**alteration**
Reported only for gene ‘labeled ‘variant’ in ‘workflow.txt’ and returns amino acid changes with corresponding positions, including replacement or absence of translational termination with ‘early_stop’ or ‘no_stop’ or ‘None’ if no variations are detected for all positions listed in ‘positions of interest’ of workflow.txt, return all variations if no positions of interest are defined or ‘ND’ if reference sequence have less than 100% coverage.

**mutation**
Reported for gene labeled ‘mutant’ in ‘workflow.txt’ and returns amino acid changes with corresponding positions or ‘None’ if no variations are detected for all positions listed in ‘positions of interest’ of workflow.txt, return ‘ND’ if no positions of interest are defined or if reference sequence have less than 100% coverage or undetected.

**modification**
Reported for gene labeled ‘regulator’ in ‘workflow.txt’ and returns nucleotide variations with corresponding positions or ‘None’ if no variations are detected for positions listed in workflow.txt, return ‘ND’ if no positions of interest are defined or ‘ND’ if reference sequence have less than 100% coverage or undetected.

**description**
Returns gene description as defined in workflow.txt

**report_type**
Returns ‘required’ or ‘conditional’ and reporting position as defined in workflow.txt

**coverage**
Reports the proportion of the reference sequence covered by a minimum acceptable depth of reads (see option –c setting before running gene_finder)

**similarity**
Reports percentage of nucleotide identities relative to the reference sequence

**Depth**
Reports the average, minimum and maximum depth of reads detected across all nucleotide positions of the gene with values separated by colons
Example value=30.59:15:55; reporting an average of 30.59, minimum of 15 and a maximum of 55 reads detected across the total mapped positions of the gene

**coverage_distribution**
Reports sections of the reference sequence covered by a minimum acceptable depth of reads
Example value=[‘1-50’, ‘100-953’]; reporting coverage with an acceptable read depth from position 1 to 50 and from 100 to 953. The absence of coverage between positions 50 and 100 could suggest a low coverage compared to the average depth, a high level of variability or deletion.

**insertions**
Reports insertion types and their corresponding positions
Example value=68:+3GGT; reporting an insertion of three nucleotides ‘GGT’ at reference sequence position 68

**deletions**
Deletion types and corresponding positions
Example value= 842:-1A
Example reporting a deletion of one nucleotide ‘A’ at reference sequence position 842

**mix**
Positions with mixed nucleotides, their natures and relative percentages
Example value=="108:C-T-0.4-0.6” ; reporting mixture of C and T with a proportion of 40 and 60%, respectively amongst a total of 120 reads mapped at position 108. 

**large_indels**
Positions with potentials large-scale indels
Example value="[(99, '19:88.89'), (100, '15:100.0')]”, reporting the proportion of reads starting or ending (88.89 and 100%) of the total number of reads (19 and 15) mapped at position 99 and 100 and thereby indicating the presence of a possible large insertion sequence between these two positions.

**mismatch**
Positions with mismatches
Example value=573:T-A;576:C-T;608:T-Y; reporting a substitution of T to A and C to T at reference sequence positions 573 and 576 whereas in position 608 the mix of C and T is reported using iupac code
 
<h3>Reference Databases</h3>
Included in the refs direcotry are example workflows for E. coli, Salmonella and Campylobacter antimicrobial resistance detection.
