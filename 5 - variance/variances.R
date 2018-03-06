library(mbtools)
library(parallel)

options(mc.cores = 6)

get_permuted_r2 <- function(counts, meta) {
    vars <- names(meta)
    permutation <- sample(vars)
    permanova <- vegan::adonis(counts ~ ., data=meta[, permutation],
                               permutations=0)
    aovs <- permanova$aov.tab[vars, "R2"]
    aovs[is.na(aovs)] <- 0
    names(aovs) <- vars
    return(aovs)
}

counts <- as.matrix(taxa_count(ps))
perms <- mclapply(1:100, function(i) get_permuted_r2(counts, meta))

