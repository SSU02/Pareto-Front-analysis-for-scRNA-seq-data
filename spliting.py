import pandas as pd
import numpy as np


df = pd.read_csv("/Users/srisruthi/Desktop/scrna/GSE131907_Lung_Cancer_cell_annotation.txt", sep = "\t")
sample = list(np.unique(df['Sample']))
all_index = list(df["Index"])
annot = {}
for i in sample:
    temp = []
    for j in all_index:
        if i in j:
            temp.append(j)
        annot[i] = temp

        
for i in list(annot.keys()):
    pd.read_csv("/Users/srisruthi/Desktop/scrna/new.csv", sep = "\t", usecols = annot[i]).to_csv(i, index = False)
    
   

        
