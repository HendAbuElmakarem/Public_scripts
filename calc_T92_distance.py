#!/usr/bin/env python3
'''
This script calculates P, Q, and theta in alignmnets for the T92 model
P= AG + TC
Q= AT+AC+TG+GC
Theta= (AG+AC+TC+TG)/2 + GC+GG+CC

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Author: Hend Abu-Elmakarem
V1_Date: 17-06-2022
V2_Date: 06-09-2022
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Output:
In this version, the script does a pairwise comparison, and produces a .tsv file with the values
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Usage: python calc_T92_parameters.py [input_alignment_file] [alignment_type]


'''
from Bio.Seq import Seq
import sys
from Bio import SeqIO
from Bio import AlignIO
import os.path
from os.path import exists
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from itertools import combinations
import math
import pandas as pd
import re


lines= []
alns = []
records = []
output = []
#read the alignment file
alignment = AlignIO.read(sys.argv[1],"fasta")


# print(alignment.get_alignment_length())
#calculate the alignment length
aln_len = alignment.get_alignment_length()
for record in alignment:
	records.append(record.id)
#loop on the alignment and pick two sequences at a time to estimate the distance between them
for pair in list(combinations(alignment, r=2)):
	m = re.findall('id*=\s*([^_]*)', str(pair))
	m =', '.join(m)
	msa = MultipleSeqAlignment(pair)
	tot_len = msa.get_alignment_length()
	substitutions = msa.substitutions
	all_sub = substitutions.select("ATCG")
	AT = all_sub['A','T'] + all_sub['T','A']
	TC = all_sub['T','C'] + all_sub['C','T']
	AG = all_sub['A','G'] + all_sub['G','A']
	AC = all_sub['A','C'] + all_sub['C','A']
	TG = all_sub['T','G'] + all_sub['G','T']
	GC = all_sub['G','C'] + all_sub['C','G']
	GG = all_sub['G','G']
	CC = all_sub['C','C']
	def nt_frequency(nt,length):
		freq = (nt/length)
		return freq

	freq_AT = nt_frequency(AT,tot_len)
	freq_TC = nt_frequency(TC,tot_len)
	freq_AG = nt_frequency(AG,tot_len)
	freq_AC = nt_frequency(AC,tot_len)
	freq_TG = nt_frequency(TG,tot_len)
	freq_GC = nt_frequency(GC,tot_len)
	freq_GG = nt_frequency(GG,tot_len)
	freq_CC = nt_frequency(CC,tot_len)

	#calculate the parameters for Tamura 1992 model
	P = freq_AG + freq_TC
	Q = freq_AT+freq_AC+freq_TG+freq_GC
	THETA = (freq_AG+freq_AC+freq_TC+freq_TG)/2 + freq_GC+freq_GG+freq_CC
	log_val_1 = 1-(1/(2*THETA*(1-THETA)))*P-Q
	log_val_2 = 1-2*Q

	#check for possible maths error and avoid caclulating log of zero
	if log_val_1 >= 0 and log_val_2 >= 0: 

		T92 = -2*THETA*(1-THETA)*math.log(1-(1/(2*THETA*(1-THETA)))*P-Q)-((1-2*THETA*(1-THETA))/2)*math.log(1-2*Q)
		df = pd.DataFrame([[m,format(T92, ".4f"),format(THETA, ".4f"),format(P, ".4f"),format(Q, ".4f"),aln_len]])
	else:

		#append the results in a dataframe
		df = pd.DataFrame([[m,"NA",format(THETA, ".4f"),format(P, ".4f"),format(Q, ".4f"),aln_len]])

	df2 = df.to_string(header=False,index=False).replace(", ","_").replace("'","").replace(" ",",")
	fob=open(sys.argv[1].split("_")[1]+'.csv','a')
	fob.write(df2)
	fob.write("\n")
	fob.close()
