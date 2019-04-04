#sanity check
library(synapser)
synLogin()
##1 compare correlation between my counts and Thanneer's counts (logCPM)
#current counts, kelsey
cc <- read.table(synGet("syn17015639")$path, sep = '\t', header=TRUE)
rownames(cc) <- cc[,1]
cc <- cc[,-1]
#former counts, thanneer
fc <- read.table(synGet("syn16847979")$path, sep = '\t', header=TRUE)
rownames(fc) <- fc[,1]
fc <- fc[,-1]

intsct <- intersect(colnames(cc), colnames(fc))
genes <- intersect(rownames(cc), rownames(fc))

c <- cc[genes,intsct]
f <- fc[genes,intsct]

compare <- cor(c,f, use = "complete.obs")
viz <- image(compare)

viz

#sanity check DLPFC 

#1.Sample to sample, and gene to gene correlations and plot the distribution of correlation coefficients in both the cases
#2. L1 or L2 Norm of the difference between your residual matrix and mine
#3. Compare the resulting covariates matrix between your and mine for quick metadata discrepancy.

