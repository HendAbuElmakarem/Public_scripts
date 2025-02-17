import os
import numpy as np
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import MutableSeq
from itertools import zip_longest

#This script concatenates multiple sequence alignments into a supermatrix, 
#it was tested on Single Copy orthologs produced using Orthofinder and aligned with MAFFT, 
#the supermatrix was then used for phylogenetic analysis. 

#list FASTA files to be concatenated
aln_files = sorted(glob.glob('OG*_new.fa'))

#Specify the species/prefix of sequences to be concatenated
raw_aln = SeqIO.parse('sps_list.fasta','fasta')

list_1 = []
files_list = []
id_seq = {}
for seq in raw_aln:
	list_1.append(str(seq.id))
	for i in list_1:
		id_seq[i]= []
for aln_file in aln_files:
	with open(aln_file,'r') as f2:
		lines = f2.readlines()
		count = 0
		for l in lines:
			if ">" in l :
				count+=1
		#add the number of sequences present in each file to be concatenate (e.g. 22 species)
		if count == 22:
			files_list.append(aln_file)
counter=0
for aln_file in aln_files:
	if aln_file in files_list:
		counter+=1
		alignment = SeqIO.parse(aln_file,'fasta')
		for rec in alignment:
			for key, val in id_seq.items():
				if key == str(rec.id.split("_")[0]):
					id_seq[key].append(str(rec.seq))

for k, v in id_seq.items():
	id_seq[k] = ''.join(id_seq[k])
with open("conc_aln.fasta","a") as f:
	for k, v in id_seq.items():
		f.write(">"+k+"\n"+v+"\n")

con_alignment = AlignIO.read("conc_aln.fasta",'fasta')
print(con_alignment.get_alignment_length())
