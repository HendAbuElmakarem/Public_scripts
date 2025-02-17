
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from itertools import zip_longest
import itertools
import glob
import os
import shutil

#This script maps protein alignment from mafft on nucleotide alignments using CDS files
#This is based on output from orthofinder, so each mafft and CDS file should have the same orthogroup number to be matched and mapped against one another 

#set the path to the working directory and the output directory, curretly set to current working directory
dir_name= os.getcwd()
final_output = os.getcwd()
os.chdir(dir_name)
#list all cds sequence files to be aligned
cds_files = sorted(glob.glob( 'cds*'))

mafft_files = sorted(glob.glob( 'mafft*'))

for (cds_file,mafft_file) in zip_longest(cds_files,mafft_files):
    #check that the sufix is the same in both files, i.e. we are dealing with the same orthogroup
    if cds_file.split('_')[0] == mafft_file.split("_")[1]:
        alignment = AlignIO.read(mafft_file,"fasta")
        alignment_len = alignment.get_alignment_length()
        pgon_cds_aln = SeqIO.parse(cds_file,"fasta")
        pgon_gene = SeqIO.parse(cds_file,"fasta")

        for cds in pgon_cds_aln:
            for record in alignment:
            #check that the species exists in both files
                if cds.id.split('_')[0] == record.id.split("_")[0]:
                    protein = []
                    gene  = []
                    #Remove stop codons from CDS files if they exist
                   if cds.seq.endswith('TGA') or cds.seq.endswith('TAG') or cds.seq.endswith('TAA'):
                       alignment_len = len(cds.seq)-3
                   else:
                       alignment_len = alignment.get_alignment_length()

                    for i in range(0,alignment_len,1):
                        protein.append(record[i])

                    for x in range(0,len(cds.seq),3):
                        gene.append(cds.seq[x:x+3])

                    for z, y in enumerate(protein):
                        if protein[z] == "-":
                            gene.insert(z,Seq("---"))
                    #The next steps ensure that we are not re-writing the newly aligned CDS files 
                    new_path = os.getcwd()
                    if os.path.isfile(os.path.join(new_path,cds_file.split('_')[1].split('.')[0]+"_cds_aln.fa")) == False:
                        with open(cds_file.split('_')[0].split('.')[0]+"_cds_aln.fa","a") as out:
                                new_seq = Seq('')
                                for n in gene:
                                        new_seq +=n
                                out.write(">"+str(record.id) + "\n" + str(new_seq)+ "\n")
for file in glob.glob("*cds_aln*"):
    shutil.move(file,new_path)

