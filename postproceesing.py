import pandas as pd

df = pd.read_csv("/Users/srisruthi/Documents/PCAfiles/EFFUSION12_PCAfile.csv",names=["x", "y", "z","a"])
#print(df)
ARC1 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/EFFUSION_12_ARC1.csv",names=["x", "y", "z"])
ARC2 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/EFFUSION_12_ARC2.csv",names=["x", "y", "z"])
ARC3 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/EFFUSION_12_ARC3.csv",names=["x", "y", "z"])
ARC4 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/EFFUSION_12_ARC4.csv",names=["x", "y", "z"])
ARC5 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/EFFUSION_12_ARC5.csv",names=["x", "y", "z"])
#ARC6 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/LUNG_T31_ARC6.csv",names=["x", "y", "z"])
#ARC7 = pd.read_csv("/Users/srisruthi/Documents/ARCfiles/LUNG_T19_ARC7.csv",names=["x", "y", "z"])
#print(df2)
mylist = [ARC1,ARC2,ARC3,ARC4,ARC5]
names = list()
answer = list()
index = 0
for j in mylist:
    print("Index:",index)
    index+=1
    names = list()
    print(len(j))
    for i in range(j.shape[0]):
        ele=j.iloc[i]
        ele_x = ele["x"]
        ele_y = ele["y"]
        ele_z = ele["z"]
        answer = (df.loc[(df["y"] == ele_x) & (df["z"] == ele_y) & (df["a"]==ele_z)])
        if (answer.empty != True):
            names.append(answer["x"].values.tolist() )


    name = [''.join(ele) for ele in names]
    print(len(name))
    new = pd.read_csv("EFFUSION_12.csv",usecols = name)
    gen = pd.read_csv("genes.csv")
    gn = pd.DataFrame(gen["Index"])
    frame = [gn,new]
    res = pd.concat(frame,axis =1)
    res.to_csv("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/P1012/EFFUSION_12/EFFUSION_12_ARC" + str(index) + ".csv", index = False)
    print("ARC " + str(index) + " DONE")
