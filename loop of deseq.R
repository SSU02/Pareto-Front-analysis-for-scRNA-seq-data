library(DESeq2)
library(BiocParallel)
param <- MulticoreParam(workers = 3, progressbar = TRUE)
register(param)
arc1<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc1.csv', header = TRUE, sep = ",", check.names = FALSE)
arc2<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc2.csv', header = TRUE, sep = ",", check.names = FALSE)
arc3<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc3.csv', header = TRUE, sep = ",", check.names = FALSE)
arc4<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc4.csv', header = TRUE, sep = ",", check.names = FALSE)
arc5<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc5.csv', header = TRUE, sep = ",", check.names = FALSE)
mylist <- list("ARC1"=arc1,"ARC2"=arc2,"ARC3"=arc3,"ARC4"=arc4,"ARC5"=arc5)
for(i in 1:5){
  for (j in (i+1):5){
    df1 <- data.frame(mylist[i])
    df1_1<- data.frame(mylist[j])
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
    dds <- DESeq(dds, fitType = "mean",parallel = TRUE, BPPARAM = param)
    res <- results(dds)
    res <- res[order(-res$padj<0.01),]
    write.csv(res,file = paste("LUNG_N06_",names(mylist[i]),"&",names(mylist[j]),"_dge.csv",sep=""))
  }
}
  