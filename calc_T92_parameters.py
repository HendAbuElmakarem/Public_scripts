#!/usr/bin/env python3
'''
This script calculates P, Q, and theta in alignmnets for the T92 model
P= AG + TC
Q= AT+AC+TG+GC
Theta= (AG+AC+TC+TG)/2 + GC+GG+CC

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
							Author: Hend Abu-Elmakarem
							Date: 17-06-2022
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Usage: python calc_T92_parameters.py [input_alignment_file] [alignment_type]


'''
import sys
from Bio import AlignIO

alignments = AlignIO.parse(sys.argv[1],sys.argv[2])

for alignment in alignments:
	tot_len = alignment.get_alignment_length()
	substitutions = alignment.substitutions
	print(substitutions)
AT = substitutions['A','T'] + substitutions['T','A']
TC = substitutions['T','C'] + substitutions['C','T']
AG = substitutions['A','G'] + substitutions['G','A']
AC = substitutions['A','C'] + substitutions['C','A']
TG = substitutions['T','G'] + substitutions['G','T']
GC = substitutions['G','C'] + substitutions['C','G']
GG = substitutions['G','G'] 
CC = substitutions['C','C']

freq_AT = AT/tot_len
freq_TC = TC/tot_len
freq_AG = AG/tot_len
freq_AC = AC/tot_len
freq_TG = TG/tot_len
freq_GC = GC/tot_len
freq_GG = GG/tot_len
freq_CC = CC/tot_len

P = freq_AG + freq_TC
Q = freq_AT+freq_AC+freq_TG+freq_GC
THETA = (freq_AG+freq_AC+freq_TC+freq_TG)/2 + freq_GC+freq_GG+freq_CC

print("P =",format(P, ".4f"))
print("Q =", format(Q, ".4f"))
print("theta =", format(THETA, ".4f"))
