library(edgeR)
raw = read.delim('/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/P0006/LUNG_N06/LUNG_N06_ARC1.csv', 
                 sep = ",", check.names=F, row.names=1)
saw = read.delim('/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/P0006/LUNG_N06/LUNG_N06_ARC2.csv', 
                 sep = ",", check.names=F, row.names=1)
rawdata = cbind(raw,saw)
df1_11 <- data.frame(colnames(raw))
df2_11 <- data.frame(colnames(saw))
meta1<- cbind(df1_11,"Archetype1")
meta2 <- cbind(df2_11,"Archetype2")
x <- c("Sample","Condition")
colnames(meta1) <- x
colnames(meta2) <- x
metaData <- rbind (meta1,meta2)
group = metaData[,-1]
dge=DGEList(rawdata, group=group)
dge=estimateDisp(dge)
dex=exactTest(dge, pair=c('Archetype1','Archetype2'))
write.csv(dex, "/Users/srisruthi/Desktop/RD lab summer intern/DE analysis/LUNG_N06_arc1&2_edgr.csv")