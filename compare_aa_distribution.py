'''
compare_aa_distribution.py

Purpose:
	Compare how site-specific amino acid frequencies differ in a particular gene between
	two user-supplied groups of samples

Dependencies
	- diamond (https://anaconda.org/bioconda/diamond)
	- mafft
	- dit
	
Usage:
	python compare_aa_distribution.py <database_path> <sample_fastq_path> <nMA_path>


Input:
	1. A protein database of the gene you want to examine (dmnd)
	2. Quality metagenomic reads (fasta)

Workflow:
	1. Creates a Master Alignment(MA) from protein database
	2. Number positions in MA (can also remove sites with low coverage)
	3. Pull out metagenomic reads mapping to database via DIAMOND
	4. From best hit, associate each position in each read with its corrosponding
		position in the numbered MA
	5. Construct an amino acid frequency table for each sample
	6. Compare amino acid frequency tables between groups

Outputs:
	- A count table showing the counts of each amino acid at each position within
		the numbered master alignment
	- A frequency table of the count table
	

-Z

'''

# Dependencies

from sys import argv
import subprocess
import dit


# Functions

def fasta_fixer(input_file,output_file):
	'''Reformat fasta file to remove newline characters within sequence

	Keyword arguments:
		input -- fasta file to fix location
		output -- fixed fasta file location
	
	'''
	out = open(output_file,'w')

	for i,l in enumerate(open(input_file,'U')):
		if l[0] == '>':
			if i == 0:
				out.write(l)
			else:
				out.write('\n'+l)
		else:
			out.write(l.strip())
	out.close()
	
def fastq_to_fasta(input_file):
	output_filename = ".".join(input_file.split(".")[:-1]) + ".fa"
	out = open(output_filename, "w")
	for i,l in enumerate(open(input_file, "U")):
		if l.startswith("@") and i % 4 == 0:
			l = l.split(" ")[0].replace("@",">")
			out.write(l + '\n')
		elif i % 4 == 1:
			out.write(l)
	out.close()
			


database_path = argv[1]
sample_path = argv[2]
sample_name = ".".join(sample_path.split(".")[:-1])
nMA_input = argv[3]


# Align all reads in the sample to MA and convert to readable SAM
# diamond_all_files.py

fastq_to_fasta(sample_path)
sample_fasta = ".".join(sample_path.split(".")[:-1]) + ".fa"

diamond_results = sample_name + '_DIAMOND_results'

subprocess.call('diamond blastx -d %s -q %s -a %s -k 1' % \
			(database_path, sample_fasta, diamond_results), shell = True)


subprocess.call('diamond view -a %s.daa -o %s.sam -f sam' % \
				(sample_name+'_DIAMOND_results', diamond_results),\
				 shell=True)


# Convert SAM file to amino acid count and frequency tables using nMA
# start from .pir file


# Grabs the length of sequences in the MA. Can probably do that when first creating the MA

align_length = 0
for l in open(nMA_input, 'U'):
	align_length = int(l.rstrip())
	break

aa_dict = {"A":0, "C":0, "D":0, "E":0, "F":0, "G":0, "H":0, "I":0, "K":0, "L":0, "M":0, "N":0, "P":0, "Q":0, "R":0, "S":0, "T":0, "V":0, "W":0, "Y":0}
aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


# Create a list the length of the MA. At each index there is a dictionary which will
# count the occurrence of each animo acid identity for all reads

count_list = []
for i in range(align_length):
	count_list.append(aa_dict.copy())

for diamond_i, diamond_l in enumerate(open(diamond_results + '.sam', 'U')):
	if diamond_i > 4: # First five lines of diamond file can be disregarded
		diamond_l = diamond_l.split('\t')
		metagenome_sequence = diamond_l[9] # index of the read sequence
		besthit_name = diamond_l[2] # name of best hit the read had in the MA
		nMA = open(nMA_input, 'U')
		for nMA_l in nMA:
			if nMA_l.startswith('>') and besthit_name in nMA_l: # find the best hit within the nMA file
				nMA_seq = next(nMA) # grab the actual sequence of the best hit
				nMA_numbering = next(nMA) # the positional numbering of the best hit
				nMA_numbering = nMA_numbering.rstrip().split(',')
				break
		tmp_sample_input_fasta_name = sample_name + '_tmp_mafft_input.fasta'
		tmp_mafft_input_file = open(tmp_sample_input_fasta_name,'w')
		tmp_mafft_input_file.write('>databaseSeq\n%s>metagenomicSeq\n%s' % \
								(nMA_seq.replace('-',''), metagenome_sequence))
		tmp_mafft_input_file.close()
		
		# Align the metagenomic read with its best hit and collects the result
		tmp_mafft_output_name = sample_name + '_unfixed_tmp_mafft_output.fasta'
		subprocess.call("mafft %s > %s" % (tmp_sample_input_fasta_name, tmp_mafft_output_name), shell = True)
		tmp_fixed_mafft_output_name = sample_name + '_tmp_mafft_output.fasta'
		fasta_fixer(tmp_mafft_output_name, tmp_fixed_mafft_output_name) # convert to single-line format
		
		# Grab both sequences from the resulting alignment
		for i, mafft_l in enumerate(open(tmp_fixed_mafft_output_name, 'U')):
			if i == 1:
				aligned_MA_seq = mafft_l # best hit sequence
			elif i == 3:
				aligned_metagenome_seq = mafft_l # sequence from the metagenome
		
		# For each spot in the alignment where neither the nMA seq or metagenomic seq
		# contained a gap, increment the amino acid identity at that position
		counter = -1
		for i, pos in enumerate(aligned_MA_seq.rstrip()):
			if pos != '-':
				counter += 1
				if aligned_metagenome_seq[i] != '-':
					identity = aligned_metagenome_seq[i] # amino acid identity at that position within the best hit
					renumbered_pos = nMA_numbering[counter] # grab associated nMA position
					dict_for_pos = count_list[int(renumbered_pos)]
					dict_for_pos[identity] += 1 # increment AA in the dictionary at the grabbed nMA position


# Scans through count_list dictionary and creates freq table
output_firstline = '\t'.join(['#POSITION', 'WT', 'SITE_ENTROPY']) + '\t' + '\t'.join(['PI_' + x for x in aa_list]) + '\n'

count_output_file = open(sample_name + '_aaCount_table.txt', 'w')
count_output_file.write(output_firstline)
freq_output_file = open(sample_name + '_aaFreq_table.txt', 'w')
freq_output_file.write(output_firstline)


for position, position_dict in enumerate(count_list):
	out_string_prefix = str(position + 1) +'\tA\t'
	count_string = ''
	for amino_acid, count in sorted(position_dict.items()):
		count_string = count_string + str(count) + '\t'
	count_string_as_list = [int(x) for x in count_string.rstrip().split('\t')]
	count_total = sum(count_string_as_list) # for making frequency table
	if count_total == 0:
		zero_line = out_string_prefix + '0' + ('\t' + str(0)) * 20 + '\n'
		count_output_file.write(zero_line)
		freq_output_file.write(zero_line)
	else:
		freq_list = [x/count_total for x in count_string_as_list]
		dit_distribution = dit.Distribution(aa_list, freq_list)
		shannon_entropy = dit.shannon.entropy(dit_distribution)
		out_firstcols = out_string_prefix + str(shannon_entropy) + '\t'
		count_output_file.write(out_firstcols + '\t'.join([str(x) for x in count_string_as_list]) + '\n')
		freq_output_file.write(out_firstcols + '\t'.join([str(x) for x in freq_list]) + '\n')


# Clean-up

subprocess.call('rm ' + tmp_sample_input_fasta_name, shell = True)
subprocess.call('rm ' + tmp_mafft_output_name, shell = True)
subprocess.call('rm ' + tmp_fixed_mafft_output_name, shell = True)

count_output_file.close()
freq_output_file.close()
		
			
					 
