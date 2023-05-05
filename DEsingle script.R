library(DEsingle)
library(BiocParallel)
param <- MulticoreParam(workers = 3, progressbar = TRUE)
raw = read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc1.csv', 
                 sep = ",", check.names=F, row.names=1)
saw = read.delim('/Users/srisruthi/Desktop/RD lab summer intern/forDGE/LUNG_N06/LUNG_N06_arc2.csv', 
                 sep = ",", check.names=F, row.names=1)
rawdata = cbind(raw,saw)
group <- factor(c(rep(1,488), rep(2,32)))
resul <- DEsingle(counts = rawdata, group = group,parallel = TRUE, BPPARAM = param)
resul.classified <- DEtype(results = resul, threshold = 0.05)
resul.sig <- resul.classified[resul.classified$pvalue.adj.FDR < 0.05, ]
write.csv(resul, "LUNG_N06_arc1&2_desingle_significant.csv")