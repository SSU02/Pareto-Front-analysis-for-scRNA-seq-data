import pandas as pd
from multiprocessing import Pool

df = pd.read_csv("/Users/srisruthi/Desktop/scrna/GSE131907_Lung_Cancer_cell_annotation.csv")

def readcsv(f):
    cell_type = {"B lymphocytes":["Index"], "Endothelial cells":["Index"], "Epithelial cells":["Index"], "Fibroblasts":["Index"], "MAST cells":["Index"], "Myeloid cells":["Index"],"T/NK cells":["Index"], "nan":["Index"]}
    for i in range(len(df)):
        if f in df["Index"][i]:
            ct = str(df["Cell_type.refined"][i])
            cell_type[ct].append(str(df["Index"][i]))
    for j in list(cell_type.keys()):
        if len(cell_type[j]) > 0:
            fname = f + "_" + j + ".csv"
            fname = fname.replace("/","_")
            fname = fname.replace(" ","_")
            pd.read_csv("/Users/srisruthi/Desktop/scrna/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt", sep = '\t', usecols=cell_type[j]).to_csv(fname, index=False)
            print(f + " DONE")

if __name__ =='__main__':
    file_name = []
    for i in df["Index"]:
        c = i.find("_")
        if i[c+1:len(i)] not in file_name:
            file_name.append(i[c+1:len(i)])
    pool = Pool(processes=8)
    pool.map(readcsv, file_name)
    
            
