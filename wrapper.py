import os #necessary dependencies
import argparse
import sys
from Bio import SeqIO

#function to parse command line arguments
def check_arg(args=None):
	parser = argparse.ArgumentParser(description='Pipeline Project Wrapper')
	parser.add_argument('-i', '--input',
		help='path to first input file',
		required='True'
		)
	parser.add_argument('-i2', '--input2',
		help='path to second input file',
		required='True'
		)
	parser.add_argument('-o', '--output',
		help='output txt file log',
		default="PipelineProject_Rumyr_Sobrepena"
		)
	return parser.parse_args(args)
#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile1 = args.input
infile2 = args.input2
outputDirectory = args.output

#making pipeline project directory to move into it
os.system(f"mkdir -p " + outputDirectory)
os.chdir(outputDirectory) 

logFile = "PipelineProject.log" #assigning txt log file 

SRR1 = infile1 #assigning first input file
SRR2 = infile2 #assigning second input file

#question 2
#Objective: Performing referenced genome mapping through bowtie2

#obtain number of reads before bowtie2 filtering
linesSRR1 = int(os.popen("wc -l < ../" + SRR1).read().strip())
linesSRR2 = int(os.popen("wc -l < ../" + SRR2).read().strip())
readsSRR1 = linesSRR1 // 4
readsSRR2 = linesSRR2 // 4
#Obtaining number of reads pairs before Bowtie2 filtering within each assembly directory 

#Retrieve the complete HCMV genome FASTA using datasets to build the index
#Using the NCBI Genome search tab, search for "Human cytomegalovirus" (HCMV) or "Human herpes virus 5"
#Obtain the command to retrieve the HCMV genome and run on the command line 
getHCMVgenome = "datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report"
os.system(getHCMVgenome)

#Unzip the dataset
unzip = "unzip ncbi_dataset.zip"
os.system(unzip)

#Use the bowtie2-build command and specify the fasta file and index name to create an index to map to for each donor
buildDatabase = "bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV"
os.system(buildDatabase)

#using bowtie to map sequences to HCMV
SRR1mapping = "bowtie2 --quiet -x HCMV -1 ../" + SRR1 + "_1.fastq -2 ../" + SRR1 + "_2.fastq -S HCMVmap1.sam"
SRR2mapping = "bowtie2 --quiet -x HCMV -1 ../" + SRR2 + "_1.fastq -2 ../" + SRR2 + "_2.fastq -S HCMVmap2.sam"
#-x refers to index we are mapping to
#-1 refers to first paired-end being mapped
#-2 refers to second paired-end being mapped
#-S refers to name of output file in .sam format
#--al-conc-gz writes paired end reads that map
os.system(SRR1mapping)
os.system(SRR2mapping)

#Count number of read pairs after Bowtie 2 filtering for each donor
afterBowtie2Filter1Reads = os.popen("samtools view -c -F 4 -f 3 HCMVmap1.sam\n").read().strip()
afterBowtie2Filter2Reads = os.popen("samtools view -c -F 4 -f 3 HCMVmap2.sam\n").read().strip() 

#converting sam files to bam
newMap1 = os.popen("samtools view -b -F 4 HCMVmap1.sam > HCMVmap1.bam").read().strip()
newMap2 = os.popen("samtools view -b -F 4 HCMVmap2.sam > HCMVmap2.bam").read().strip()

#converting bam files to fastq
SRR1Mapped = os.popen("samtools fastq -f 3 -1 " + SRR1 + "_mapped_1.fq.gz -2 " + SRR1 + "_mapped_2.fq.gz HCMVmap1.bam").read().strip()
SRR2Mapped = os.popen("samtools fastq -f 3 -1 " + SRR2 + "_mapped_1.fq.gz -2 " + SRR2 + "_mapped_2.fq.gz HCMVmap2.bam").read().strip()

#making and editing output file
with open(logFile, 'a') as f:
	f.write("Donor 1 (2dpi) had " + str(readsSRR1) + " read pairs before Bowtie 2 filtering and " + afterBowtie2Filter1Reads + " read pairs after.\n")
	f.write("Donor 2 (6dpi) had " + str(readsSRR2) + " read pairs before Bowtie 2 filtering and " + afterBowtie2Filter2Reads + " read pairs after.\n\n")

#question 3
#Objective: Using the Bowtie2 output reads, assemble both transcriptomes together to produce assembly via SPAdes. You should use a k-mer size of 99.
#including the command in the pipeline log
spadesCommand = os.popen("spades.py -k 99 -t 2 --only-assembler -1 ./" + SRR1 + "_mapped_1.fq.gz -2 ./" + SRR1 + "_mapped_2.fq.gz -1 ./" + SRR2 + "_mapped_1.fq.gz -2 ./" + SRR2 + "_mapped_2.fq.gz -o " + SRR1 + "_and_" + SRR2 + "_assembly/").read().strip()
spadesC = f"spades.py -k 99 -t 2 --only-assembler -1 ./" + SRR1 + "_mapped_1.fq.gz -2 ./" + SRR1 + "_mapped_2.fq.gz -1 ./" + SRR2 + "_mapped_1.fq.gz -2 ./" + SRR2 + "_mapped_2.fq.gz -o " + SRR1 + "_and_" + SRR2 + "_assembly/\n\n"
with open(logFile, 'a') as f: #writing spades command to log file
	f.write(spadesC)
os.system(spadesCommand)
	
#question 4
contigs = [] #making a list of contig lengths
contigs1000 = 0
bpAssembly = 0

#getting number of contigs over 1000 bp and length of bpAssembly
for record in SeqIO.parse(SRR1 + "_and_" + SRR2 + "_assembly/contigs.fasta", "fasta"):
	contigLength = len(record.seq)
	contigs.append((record.id, contigLength, str(record.seq)))
	if contigLength > 1000:
		contigs1000 += 1
		bpAssembly += contigLength

#writing output to pipeline log
with open(logFile, 'a') as f:
	f.write("There are " + str(contigs1000) + " contigs > 1000 bp in the assembly \n")
	f.write("There are " + str(bpAssembly) + " bp in the assembly.\n\n")

#Question 5
#retrieve the longest contig from SPAdes assembly
#make a local database from the Betaherpesvirinae subfamily
#Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily; should only keep the best alignment (HSP) for any single query-subject pair of sequences
#For the top 10 hits, write to log file: Subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit scoree, e-value, and subject title

#Obtain the Betaherpesvirinae sequences from NCBI using the command line
getBHVgenome = "datasets download genome taxon 'Betaherpesvirinae' --reference" 
os.system(getBHVgenome)
os.system(unzip) #unzip current ncbi_dataset.zip

#combine all betaherpesvirinae sequences into single file 
catFNA = "cat ncbi_dataset/data/*/*genomic.fna > betaherpesvirinae.fna"
os.system(catFNA)

#Making local database from the Betaherpesvirinae subfamily
makeLocalDB = "makeblastdb -in betaherpesvirinae.fna -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl"
#-in specifies input file
#-out specifies output file
#-title specifies desired title
#-dbtype specifies type of database
os.system(makeLocalDB)

#getting longest contig from contig list
contigs.sort(key=lambda x: x[1], reverse=True)
longestContigID, longestContigLength, longestContigSequence = contigs[0]

#writing longest contig into fasta file readable for blast
with open("longestContig.fasta", 'w') as f: 
	f.write(f">" + str(longestContigID) + "\n" + str(longestContigSequence) + "\n\n")

#running blastn to compare longest contig with betaherpesvirinae sequences
blastCMD = 'blastn -query longestContig.fasta -db betaherpesvirinae -out betaherpesvirinae_blastn_results.tsv -max_target_seqs 10 -max_hsps 1 -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
#outfmt describes format of the output
#6 refers to tabular
#sacc refers to subject accession
#pident refers to percent identity
#length refers to alignment length
#qstart refers to start of alignment in query
#qend refers to end of alignment in query
#sstart refers to start of alignment in subject
#send refers to end of alignment in subject
#bitscore refers to bit score
#evalue refers to expected value
#stitle refers to subject title 
#max_target_seqs describes most sequences wanted
#max_hsps refers to keeping best alignment
os.system(blastCMD)

with open(logFile, 'a') as f: #writes results of blast to Pipeline log
	f.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
	with open("betaherpesvirinae_blastn_results.tsv", 'r') as g:
		f.write(g.read())




