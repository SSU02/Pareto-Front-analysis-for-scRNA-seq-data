library(DEsingle)
library(BiocParallel)
param <- MulticoreParam(workers = 3, progressbar = TRUE)
register(param)
arc1= read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc1.csv', sep = ",", check.names=F, row.names=1)
arc2= read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc2.csv', sep = ",", check.names=F, row.names=1)
arc3= read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc3.csv', sep = ",", check.names=F, row.names=1)
arc4= read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc4.csv', sep = ",", check.names=F, row.names=1)
arc5= read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc5.csv', sep = ",", check.names=F, row.names=1)
mylist <- list("ARC1"=arc1,"ARC2"=arc2,"ARC3"=arc3,"ARC4"=arc4,"ARC5"=arc5)
for(i in 1:5){
  for (j in (i+1):5){
    raw <- data.frame(mylist[i])
    saw <- data.frame(mylist[j])
    rawdata = cbind(raw,saw)
    group <- factor(c(rep(1,length(data.frame(mylist[i]))), rep(2,length(data.frame(mylist[j])))))
    resul <- DEsingle(counts = data.frame(rawdata), group = group, parallel = TRUE, BPPARAM = param)
    resul.classified <- DEtype(results = resul, threshold = 0.05)
    resul.sig <- resul.classified[resul.classified$pvalue.adj.FDR < 0.05, ]
    write.csv(resul.sig, file = paste("LUNG_N06_",names(mylist[i]),"&",names(mylist[j]),"_DEsingle_signif.csv",sep=""))
 }
}  