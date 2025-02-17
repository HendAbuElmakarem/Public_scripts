import itertools
import re
import glob
import os

#This scripts extracts the values of dN, dS and omega of YN00 and LPB93 from PAML output
#The script uses regex to find the desired values in the text  
dn1="dN +"
result_files = glob.glob("yn_*txt")

with open("dN_dS_pairs_output.txt","w+") as f11:
     f11.write("Method"+ "\t" + "Orthogroup" + "\t" + "Pair" + "\t" + "dN" + "\t" + "dS" + "\n")

for result_file in result_files:
    f = open(result_file) 
    iter_file = iter(f)
    with open("dN_dS_pairs_output.txt","a+") as f1:
         for line in f:
             if dn1 in line:
                line = next(iter_file)
                line = next(iter_file)
                parts = line.split()
                if len(parts)> 6:
                   f1.write("YN00" + "\t" + result_file.split("_")[3].replace(".txt"," ") + "\t" + result_file.split("_")[1] + "_" + result_file.split("_")[2] + "\t" + parts[7] + "\t" + parts[10] + "\n")
             elif line.startswith("LPB93"):
                parts = line.split()
                f1.write("LPB93" + "\t" + result_file.split("_")[3].replace(".txt"," ") + "\t" + result_file.split("_")[1] + "_" + result_file.split("_")[2] + "\t" + parts[6] + "\t" + parts[3] + "\n")
