# myrepo 
# Objective: Retrieving sequence data from NCBI
 Step 1: Go to NCBI link for Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360 \n
 Step 2: Click the link to the Run ID that corresponds to the transcriptome
 Step 3: Click the Data access tab
 Step 4: Copy the link of the SRA Normalized data using the command line with the wget command, which downloads the file associated with the link provided
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030


# Step 5: Go to NCBI link for Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov
# Repeat steps 2-4 
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033 

# Objective: Convert transcriptomes to paired-end fastq files 

# Step 1: To uncompress the data, use faster-qdump from NCBI's SRA-Toolkit
fasterq-dump SRR5660030
fasterq-dump SRR5660033
# fasterq-dump converts SRA file to FASTQ format and results in two unpaired reads files being written into two fastq files

# Step 2: To create paired end reads of each data
spades.py -k 77,99,127 -t 2 --only-assembler -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -o SRR5660030_assembly/
spades.py -k 77,99,127 -t 2 --only-assembler -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -o SRR5660033_assembly/ 
# -k is list of k-mer sizes that SPAdes uses to build de Bruijn graphs and tests which k finds the best assembly 
# -t is number of threads
# --only-assembler runs only assembling (without read error correction)
# -1 specifies first input file
# -2 specifies second inpute file 
# -o specifies the directory to store all resulting files

# output is also written to a file in the output directory called spades.log 


# Performing de novo genome assembly through spades.py


# QUESTION 2
# Objective: Performing referenced genome mapping through bowtie2

# Obtaining number of reads pairs before Bowtie2 filtering within each assembly directory 
grep -c "^>" scaffolds.fasta #do this in each assembly directory of the two donors

#
nano question1Output.py 

# Step 1: Retrieve the complete HCMV genome FASTA using datasets to build the index
# Using the NCBI Genome search tab, search for "Human cytomegalovirus" (HCMV) or "Human herpes virus 5"
# Obtain the command to retrieve the HCMV genome and run on the command line 
datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report

# Step 2: Unzip the dataset
unzip ncbi_dataset.zip

# Step 3 Use the bowtie2-build command and specify the fasta file and index name to create an index to map to for each donor
bowtie2-build ncbi_dataset/data/GCF000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV

bowtie2 -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMVmap.sam --al-conc-gz SRR5660030_mapped_%.fq.gz
#-x refers to index we are mapping to
#-1 refers to first paired-end being mapped
#-2 refers to second paired-end being mapped
#-S refers to name of output file in .sam format
#--al-conc-gz writes paired end reads that map

# Count number of read pairs after Bowtie 2 filtering for donor 1
grep -c "^@" SRR5660030_mapped_%.fq.gz
# -c counts matching lines

bowtie2 -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMVmap.sam --al-conc-gz SRR5660033_mapped_%.fq.gz
# -x refers to index we are mapping to
# -1 refers to first paired-end being mapped
# -2 refers to second paired-end being mapped
# -S refers to name of output file in .sam format
# --al-conc-gz writes paired end reads that map

# Count number of read pairs after Bowtie 2 filtering for donor 2
grep -c "^@" SRR5660033_mapped_%.fq.gz
# -c counts matching lines


# Question 3
# Objective: Using the Bowtie2 output reads, assemble both transcriptomes together to produce assembly via SPAdes. You should use a k-mer size of 99.
spades.py -k 99 -t 2 --only-assembler -1 SRR5660030_mapped_%.fq.gz -2 SRR5660033_mapped_%.fq.gz -o SRR5660030_and_SRR5660033_assembly/

# Write SPAdes command to log file
# ...

# Question 5
# Python, retrieve the longest contig from SPAdes assembly
# make a local database from the Betaherpesvirinae subfamily
# Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily; should only keep the best alignment (HSP) for any single query-subject pair of sequences
# For the top 10 hits, write to log file: Subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit scoree, e-value, and subject title

# Building off of question 4 code, which obtained longest contig 

# Step 1: Obtain the Betaherpesvirinae sequences from NCBI using the command line
datasets download virus genome betaherpesvirinae --refseq --include genome 

# Step 2: Making local database from the Betaherpesvirinae subfamily
makeblastdb -in -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl
# -in specifies input file
# -out specifies output file
# -title specifies desired title
# -dbtype specifies type of database
