#!/usr/bin/env python
# coding: utf-8

# In[15]:







import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--infile',help='The multiple sequence alignment must be in clustal format')

parser.add_argument('-o','--outfile',help='Name the output phylip file')
args = parser.parse_args()
    
infile = args.infile
outfile = args.outfile




def filter_converter(infile,outfile):
    import pandas as pd
    import phylopandas as ph
    df1 = ph.read_clustal(infile)
    
    result_df = df1.drop_duplicates(subset=['sequence'], keep='first')
    result_df.to_phylip(outfile)
if __name__ == '__main__':
    print(filter_converter(infile,outfile))


# In[ ]:




