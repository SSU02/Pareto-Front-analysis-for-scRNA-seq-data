library(corrplot)

data <- read.csv("correlation.csv", header = T, row.names = 1, sep = ",", check.names = FALSE)

res1 <- cor(data, method = "pearson")
write.csv(res1, "Pearson-coeff.csv")


res <- cor(data, method = "spearman")
write.csv(res, "spearman.csv")

pdf("Spearman.pdf")
corrplot(res, method = "pie", type = "lower")
dev.off()

pdf("Pearson.pdf")
corrplot(res1, method = "pie", type = "lower")
dev.off()           