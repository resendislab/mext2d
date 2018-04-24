library(mbtools)
library(parallel)

options(mc.cores = 6)

counts <- as.matrix(taxa_count(ps))
meta <- meta[rownames(counts), ]
D <- vegan::vegdist(counts, method="bray")
	perms <- vegan::adonis(D ~ ., data=meta, permutations=10000, parallel=TRUE)

