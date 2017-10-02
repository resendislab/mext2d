# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Scripts for running association tests on the genus level

devtools::load_all("mbtools")
library(IHW)
library(magrittr)
library(pbapply)

ps <- readRDS("data/taxonomy.rds")
counts <- as.matrix(taxa_count(ps))

# remove duplicates
s <- tstrsplit(rownames(counts), "_")[[1]]  # sample ids
dupes <- s %in% (which(table(s) != 1) %>% names()) %>% which()
counts <- counts[-dupes, ]
s <- s[-dupes]

# get meta data and remove missing samples
meta <- readxl::read_excel("data/clinical.xlsx", sheet = "clean + english")
type <- readxl::read_excel("data/clinical.xlsx", sheet = "variables")
meta <- meta %>% types(type) %>% standardize() %>% as.data.table()
setkey(meta, id)
non_missing <- which(s %in% meta$id)
s <- s[non_missing]
counts <- counts[non_missing, ]
meta <- meta[s]
meta$status <- as.integer(meta$status)

# Assemble DESeq2 data set
library(DESeq2)

logger <- file("testing.log", open="w")
sink(file = logger, type = "message")
confounders <- c("gender")
dds <- DESeqDataSetFromMatrix(t(counts), meta, design = ~ gender + status)
dds <- estimateSizeFactors(dds, type = "poscount")
dds <- DESeq(dds, parallel = TRUE)
exclude <- c("status", "id", "stool_dna", confounders)
vars <- names(meta[, !exclude, with=F])
multi <- results(dds, name = "status")
multi <- lfcShrink(dds, coef = 3, res = multi)
multi <- as.data.table(multi)
multi$genus <- colnames(counts)
multi$variable <- "status"
multi$n_test <- nrow(counts)

tests <- pblapply(vars, function(v) {
    message(paste("Testing", v, "...\n"))
    good <- !is.na(meta[[v]])
    dds <- DESeqDataSetFromMatrix(t(counts[good, ]), as.data.frame(meta[good]),
                                  design = reformulate(c(confounders, v)))
    dds <- estimateSizeFactors(dds, type = "poscount")
    dds <- DESeq(dds, parallel = TRUE)
    res <- lfcShrink(dds, coef = length(resultsNames(dds)), res = results(dds))
    res <- as.data.table(res)
    res$genus <- colnames(counts)
    res$variable <- v
    res$n_test <- sum(good)
    res
})
sink()
close(logger)

# Remove genus with a baseMean smaller 1 to avoid a bimodal pval distribution
multi <- rbind(multi, rbindlist(tests))[baseMean >= 1]
weighting <- ihw(pvalue ~ baseMean, multi)
multi[, padj := adj_pvalues(weighting)]
fwrite(multi[order(padj, variable)], "data/tests.csv")
