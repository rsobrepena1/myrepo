# myrepo 

Step 0: Ensure that the required dependencies are installed (Bowtie2, SPAdes, BLAST+, Biopython). You can use pip or conda for installation. 

#QUESTION 1  

Objective: Retrieving sequence data from NCBI  

Step 1: Go to NCBI link for Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360  
Step 2: Click the link to the Run ID that corresponds to the transcriptome  
Step 3: Click the Data access tab  
Step 4: Copy the link of the SRA Normalized data using the command line with the wget command, which downloads the file associated with the link provided  

Step 5: In the command line:  
```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030  
```

Step 6: Go to NCBI link for Donor 2 (6dpi): https://www.ncbi.nlm.nih.gov
Repeat steps 2-4 for Donor 2   
```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033 
```
if file is large, Subsample data to write first 10,000 reads to a file  
```
head -n 40000 SRR5660030 > sampleSRR5660030
head -n 40000 SRR5660033 > sampleSRR5660033
```
#Objective: Convert transcriptomes to paired-end fastq files  

Step 1: To uncompress the data, use faster-qdump from NCBI's SRA-Toolkit
fasterq-dump converts SRA file to FASTQ format and results in two read files being written into two fastq files  
```
fasterq-dump sampleSRR5660030 
fasterq-dump sampleSRR5660033
```

Paired end reads of sample data have been provided in the myrepo. 
These should be downloaded to your directory.


Complete questions 2-5 by calling wrapper.py through the terminal with 
```
python wrapper.py --input [your first file] --input2 [your second file]""
```
```
Example: python wrapper.py --input sampleSRR5660030 --input2 sampleSRR5660033 
```

Output is written to a file in the output directory: Pipeline_log