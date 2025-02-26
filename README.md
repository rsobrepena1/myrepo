# myrepo 

Step 0: Ensure that the required dependencies are installed (Bowtie2, SPAdes, BLAST+, Biopython)

#QUESTION 1

#Objective: Retrieving sequence data from NCBI

Step 1.1: Go to NCBI link for Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Step 1.2: Click the link to the Run ID that corresponds to the transcriptome
Step 1.3: Click the Data access tab
Step 1.4: Copy the link of the SRA Normalized data using the command line with the wget command

In the command line: wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030 #Donor 1 (2dpi)

wget downloads the file associated with the link provided

Step 1.5: Go to NCBI link for Donor 2 (6dpi): https://www.ncbi.nlm.nih.gov
Repeat steps 1.2-1.4 
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033 #Donor 2 (6dpi)

if file is large, Subsample data to write first 10,000 reads to a file as a sample test file
head -n 40000 SRR5660030 > sampleSRR5660030
head -n 40000 SRR5660033 > sampleSRR5660033

#Objective: Convert transcriptomes to paired-end fastq files 

Step 1: To uncompress the data, use faster-qdump from NCBI's SRA-Toolkit
fasterq-dump sampleSRR5660030 
fasterq-dump sampleSRR5660033
fasterq-dump converts SRA file to FASTQ format and results in two unpaired reads files being written into two fastq files

Step 2: To create paired end reads of each data
spades.py -k 77,99,127 -t 2 --only-assembler -1 sampleSRR5660030_1.fastq -2 sampleSRR5660030_2.fastq -o SRR5660030_assembly/ 
spades.py -k 77,99,127 -t 2 --only-assembler -1 sampleSRR5660033_1.fastq -2 sampleSRR5660033_2.fastq -o SRR5660033_assembly/
spades.py 
-k is list of k-mer sizes that SPAdes uses to build de Bruijn graphs and tests which k finds the best assembly 
-t is number of threads
--only-assembler runs only assembling (without read error correction)
-1 specifies first input file
-2 specifies second input file 
-o specifies the directory to store all resulting files

output is also written to a file in the output directory: Pipeline_log

Complete questions 2-5 by calling wrapper.py with "python wrapper.py --input [your first file] --input2 [your second file]""
Example: python wrapper.py --input SRR5660030 --input2 SRR5660033 
