library(MAST)
#reading the first file
df1 <- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc1.csv', header = TRUE, sep = ",", check.names = FALSE)
#reading the second file
df1_1<- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc2.csv', header = TRUE, sep = ",", check.names = FALSE)
#removing index column from 2nd file
df2 <- df1_1[,-1]
#binding file1 and file2 columnwise
df <- cbind(df1,df2)
#making index column names in the df file as row names
rownames(df) <- df[,1]
#removing index column from df to get df3
df3<- df[,-1]
#getting header names from file1 as a new column in df1_11
df1_11 <- data.frame(colnames(df1[,-1]))
#getting header names from file2 as a new column in df2_11
df2_11 <- data.frame(colnames(df1_1[,-1]))
#creating a meta file to give as phenotype data 
meta1<- cbind(df1_11,"Archetype1")
meta2 <- cbind(df2_11,"Archetype2")
x <- c("wellKey","condition")
colnames(meta1) <- x
colnames(meta2) <- x
metaData <- rbind (meta1,meta2)
rownames(metaData) <- metaData[,1]
#phenotype data as pd
#pd <- new("AnnotatedDataFrame", data = metaData)
#feature data as fd (getting gene names)
fdtry1 <- df1[,1]
fdtry2 <- cbind(fdtry1,"dummy","dummy1")
colnames(fdtry2) <- c("primerid","entrez","symbolid")
fdtry2 <- data.frame(fdtry2)
rownames(fdtry2) <- fdtry2[,1]
#fd <- new("AnnotatedDataFrame", data = fdtry2)
#dds <- FromMatrix(exprsArray = as.matrix(df3), cData = pd,fData = fd)
sca <- FromMatrix(exprsArray = as.matrix(df3), cData = metaData,fData = fdtry2)
cond<-factor(colData(sca)$condition)
cond<-relevel(cond,"Archetype2")
colData(sca)$condition<-cond
zlmCond <- zlm(~condition, sca)
summaryCond <- summary(zlmCond, doLRT='conditionArchetype1') 
summaryDt <- summaryCond$datatable
dummy <- filter(summaryDt, component == "logFC")
write.csv(dummy, "dummy.csv")


                  
                  
                  
                  