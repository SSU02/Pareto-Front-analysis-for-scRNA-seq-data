library(DESeq2)
library(ggplot2)
library(glmpca)

df1 <- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc1.csv', header = TRUE, sep = ",", check.names = FALSE)
df1_1<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc2.csv', header = TRUE, sep = ",", check.names = FALSE)
df2 <- df1_1[,-1]
df <- cbind(df1,df2)
df1_11 <- data.frame(colnames(df1[,-1]))
df2_11 <- data.frame(colnames(df1_1[,-1]))
meta1<- cbind(df1_11,"Archetype1")
meta2 <- cbind(df2_11,"Archetype2")
x <- c("Sample","Condition")
colnames(meta1) <- x
colnames(meta2) <- x
metaData <- rbind (meta1,meta2)
count<- round(df[,-1],0)
countData <-cbind(df[,1],count)
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Condition, tidy = TRUE)
dds <- DESeq(dds, fitType = "mean")
res <- results(dds)
summary(res)
res <- res[order(-res$padj<0.01),]
write.csv(res, "LUNG_N06_arc1&2_dge.csv")