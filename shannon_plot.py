#!/usr/bin/env python
# coding: utf-8

# In[10]:

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a','--alnfile',help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
parser.add_argument('-f','--fileformat',help='Specify the format of the input MSA to be passed in to AlignIO.')
parser.add_argument('-t','--titlename',help='The title name for the Shannon entropy plot')
parser.add_argument('-o','--outputfig',help='Name the output png file')
args = parser.parse_args()
    
alnfile = args.alnfile
fileformat = args.fileformat
titlename = args.titlename
outputfig = args.outputfig


def shannon_plot(alnfile,fileformat,titlename,outputfig):
    import pandas as pd 
    from Bio import AlignIO
    import math
    from scipy.stats import entropy
    import numpy as np
    import matplotlib.pyplot as plt
    

    alignment = AlignIO.read(alnfile, fileformat)

    for record in alignment:
        len_aln = len(record)
    
    column_list = []

#extract each column as strings
    for i in range(len_aln):
        column_list.append(alignment[:,i])
    lis = [] #list of columns
    for t in range(len(column_list)):
    
    
        seq = column_list[t].split(",")
        l = seq[0]
        splitted = []
        s = [i for i in l]
        splitted.append(s)
        lis.append(splitted)
    v = len(lis)
    site = []
    h_list = []

    for y in range(v):
        site.append(y)
        unique_base = np.unique(lis[y][0])
    
        M = len(lis[0][0])
        entropy_list = []
        for base in unique_base:
        
            n_i = lis[y][0].count(base) # Number of residues of type i
            P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
            entropy_i = P_i*(math.log(P_i,2))
            entropy_list.append(entropy_i)
        

    
        sh_entropy = -(sum(entropy_list))
        if sh_entropy == 0:
            h_list.append(0)
        else:
            h_list.append(sh_entropy)
    
 
    x = range(len(h_list))
    x_values = [i + 1 for i in x]
    plt.plot(x_values, h_list)
    plt.title(titlename)
    plt.xlabel('Position')
    plt.ylabel('Entropy')
    plt.savefig(outputfig)
    


       
if __name__ == '__main__':
    print(shannon_plot(alnfile,fileformat,titlename,outputfig))