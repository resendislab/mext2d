# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Scripts for running association tests on the genus level

library(mbtools)
library(magrittr)
library(DESeq2)

zero_fraction <- 0.5
min_abundance <- 8
skip_all <- file.exists("../data/tests.csv")

who <- function(glucose_0, glucose_120) {
    class <- rep("NGT", length(glucose_0))
    class[is.na(glucose_0) | is.na(glucose_120)] <- NA
    class[glucose_0 < 126 & glucose_120 >= 140 & glucose_120 < 200] <- "IGT"
    class[glucose_0 >= 126 & glucose_120 >= 200] <- "T2D"
    return(factor(class, level=c("NGT", "IGT", "T2D")))
}

ps <- readRDS("../data/taxonomy.rds")
counts <- as.matrix(taxa_count(ps))

# remove duplicates
s <- tstrsplit(rownames(counts), "_")[[1]]  # sample ids
dupes <- s %in% (which(table(s) != 1) %>% names()) %>% which()
counts <- counts[-dupes, ]
if (length(dupes) > 0) cat("found duplicates:", unique(s[dupes]), "\n")
s <- s[-dupes]

# get meta data and remove missing samples
meta <- readxl::read_excel("../data/clinical.xlsx", sheet = "clean + english")
meta$who <- who(meta$glucose_0, meta$glucose_120)
type <- readxl::read_excel("../data/clinical.xlsx", sheet = "variables")
meta <- meta %>% types(type) %>% standardize() %>% as.data.table()
setkey(meta, id)
non_missing <- which(s %in% meta$id)
s <- s[non_missing]
counts <- counts[non_missing, ]
rownames(counts) <- s

# Exclude taxa that are absent in more than zero_fraction samples
fraction_exclude <- colSums(counts >= 1) / nrow(counts) < zero_fraction
cat("removed genera due to missing reads:",
    sum(fraction_exclude), "\n")
counts <- counts[, !fraction_exclude]
meta <- meta[s]
meta$status <- as.integer(meta$status)

# Exclude samples added by error (status = 7)
meta <- meta[status < 6]
counts <- counts[meta$id, ]

# Assemble DESeq2 data set
confounders <- c("gender")
exclude <- c("id", "stool_dna", confounders)
vars <- names(meta[, !exclude, with = F])
tests <- NULL

good <- !is.na(meta$who)
dds <- DESeqDataSetFromMatrix(t(counts[good, ]), as.data.frame(meta[good]),
                              design = ~ gender + who)
dds <- estimateSizeFactors(dds, type = "poscount")
dds <- DESeq(dds, parallel = TRUE)

multi <- NULL
comps <- matrix(c("T2D", "T2D", "IGT", "NGT", "IGT", "NGT"), ncol=2)
for (i in 1:nrow(comps)) {
    name <- paste0(comps[i, 1], "_vs_", comps[i, 2])
    contrast <- c("who", comps[i, 1], comps[i, 2])
    res <- results(dds, contrast = contrast)
    res <- lfcShrink(dds, contrast = contrast, res = res)
    res <- as.data.table(res)
    res$genus <- colnames(counts)
    res$variable <- name
    res$n_test <- sum(good)
    multi <- rbind(multi, res)
}
# Remove genus with a small baseMean to avoid a bimodal pval distribution
multi <- multi[baseMean >= min_abundance]
multi[, padj := p.adjust(pvalue)]
fwrite(multi[order(padj, variable)], "../data/tests_who_genus.csv")
