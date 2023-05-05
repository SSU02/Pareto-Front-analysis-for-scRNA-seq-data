library(MAST)
library(dplyr)
library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
patient = list("P0034","P1006","P1011","P1012", "P1013");
P=patient[5]
s = list("EBUS_13", "EFFUSION_13", "BRONCHO_11");
sample = s[2];
arc1 <- read.csv(paste0("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/",P,"/",sample,"/",sample,"_ARC1.csv"), header = TRUE, sep = ",", check.names = FALSE)
arc2 <- read.csv(paste0("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/",P,"/",sample,"/",sample,"_ARC2.csv"), header = TRUE, sep = ",", check.names = FALSE)
arc3 <- read.csv(paste0("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/",P,"/",sample,"/",sample,"_ARC3.csv"), header = TRUE, sep = ",", check.names = FALSE)
arc4 <- read.csv(paste0("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/",P,"/",sample,"/",sample,"_ARC4.csv"), header = TRUE, sep = ",", check.names = FALSE)
#arc5 <- read.csv(paste0("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/",P,"/",sample,"/",sample,"_ARC5.csv"), header = TRUE, sep = ",", check.names = FALSE)
#arc6 <- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/P0031/LUNG_T31/LUNG_T31_ARC6.csv', header = TRUE, sep = ",", check.names = FALSE)
#arc7 <- read.csv('/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/P0019/LUNG_T19/LUNG_T19_ARC7.csv', header = TRUE, sep = ",", check.names = FALSE)
#mylist <- list("ARC1"=arc1,"ARC2"=arc2,"ARC3"=arc3,"ARC4"=arc4,"ARC5"=arc5,"ARC6"=arc6,"ARC7"=arc7)
#mylist <- list("ARC1"=arc1,"ARC2"=arc2,"ARC3"=arc3,"ARC4"=arc4,"ARC5"=arc5,"ARC6"=arc6)
#mylist <- list("ARC1"=arc1,"ARC2"=arc2,"ARC3"=arc3,"ARC4"=arc4,"ARC5"=arc5)
mylist <- list("ARC1"=arc1,"ARC2"=arc2,"ARC3"=arc3,"ARC4"=arc4)

foreach(i = 1:4) %do% {
  foreach (j = (i+1):4) %do% {
    df1 <- data.frame(mylist[i])
    df1_1<- data.frame(mylist[j])
    df2 <- df1_1[,-1]
    df <- cbind(df1,df2)
    rownames(df) <- df[,1]
    df3<- df[,-1]
    df1_11 <- data.frame(colnames(df1[,-1]))
    df2_11 <- data.frame(colnames(df1_1[,-1]))
    meta1<- cbind(df1_11,"Archetype1")
    meta2 <- cbind(df2_11,"Archetype2")
    x <- c("wellKey","condition")
    colnames(meta1) <- x
    colnames(meta2) <- x
    metaData <- rbind (meta1,meta2)
    rownames(metaData) <- metaData[,1]
    fd1 <- df1[,1]
    fd2 <- cbind(fd1,"dummy","dummy1")
    colnames(fd2) <- c("primerid","entrez","symbolid")
    fd2 <- data.frame(fd2)
    rownames(fd2) <- fd2[,1]
    sca <- FromMatrix(exprsArray = as.matrix(df3), cData = metaData,fData = fd2)
    cond<-factor(colData(sca)$condition)
    cond<-relevel(cond,"Archetype2")
    colData(sca)$condition<-cond
    zlmCond <- zlm(~condition, sca,parallel = TRUE)
    summaryCond <- summary(zlmCond, doLRT='conditionArchetype1') 
    summaryDt <- summaryCond$datatable
    res <- filter(summaryDt, component == "logFC")
    write.csv(res,file = paste0("/Users/srisruthi/Desktop/RD lab summer intern/cancer samples/",P,"/",sample,"/",sample,"_",names(mylist[i]),"&",names(mylist[j]),"_dge.csv",sep=""))
  }
}