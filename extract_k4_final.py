from Bio.Seq import Seq
import sys
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
from itertools import groupby
from Bio.Seq import MutableSeq
from itertools import zip_longest

#This script extracts fourfold degenerate sites (K4) from a multiple sequence alignment
#K4 sites are expected to be neutral as any nucleotide change at these sites does not change the resulting amino acid.

#Read the alignment
alignment = AlignIO.read(sys.argv[1],"fasta")
align_2 = AlignIO.read(sys.argv[1],"fasta")

#Specify codon at which the third position is a K4 site 
k4_codons = ['CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC',
				'GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG']
#Store possible values of gaps or Ns to be discarded from the analysis
all_pos_gap_1 = set(['---','nnn','NNN'])

alignment_len = alignment.get_alignment_length()
K4_id_matched_cols = set()

#Loop on each codon in the alignment (a step of three)
for i in range(0, alignment_len, 3):
	new_aln = alignment[:,i:i + 3]
	col = []
	for rec in new_aln:
		col.append(str(rec.seq))
	for j in range(0, len(col)):
		all_to_ignore_1 = any(elem in all_pos_gap_1 for elem in col)		
		K4_result =  all(elem in k4_codons for elem in col)

		if all_to_ignore_1:
			continue;

		#If site is K4 across the alignment then append its index to a set
		if K4_result:

			K4_id_matched_cols.add(i+2)

first_index_list = sorted(K4_id_matched_cols)
for record1, record2 in zip_longest(alignment,align_2) :
    record1.seq = MutableSeq(record1.seq)
    record2.seq = MutableSeq('')
    for item in first_index_list:
        record2.seq.append(record1.seq[item])
print(align_2.get_alignment_length())
AlignIO.write(align_2, "k4_"+ sys.argv[1], "fasta")

