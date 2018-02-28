# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Scripts for running association tests on the genus level

library(mbtools)
library(magrittr)
library(DESeq2)

zero_fraction <- 0.1
min_abundance <- 10
skip_all <- file.exists("../data/tests_family.csv")

ps <- readRDS("../data/taxonomy.rds")
counts <- as.matrix(taxa_count(ps, "family"))

# remove duplicates
s <- tstrsplit(rownames(counts), "_")[[1]]  # sample ids
dupes <- s %in% (which(table(s) != 1) %>% names()) %>% which()
counts <- counts[-dupes, ]
if (length(dupes) > 0) cat("found duplicates:", unique(s[dupes]), "\n")
s <- s[-dupes]

# get meta data and remove missing samples
meta <- readxl::read_excel("../data/clinical.xlsx", sheet = "clean + english")
type <- readxl::read_excel("../data/clinical.xlsx", sheet = "variables")
meta <- meta %>% types(type) %>% standardize() %>% as.data.table()
setkey(meta, id)
non_missing <- which(s %in% meta$id)
s <- s[non_missing]
counts <- counts[non_missing, ]
rownames(counts) <- s

# Exclude taxa that are absent in more than zero_fraction samples
fraction_exclude <- (colSums(counts >= 1) / nrow(counts)) < zero_fraction
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


if (!skip_all) {
    for (v in vars) {
        cat("---\nTesting", v, "\n---\n")
        good <- !is.na(meta[[v]])
        dds <- DESeqDataSetFromMatrix(t(counts[good, ]),
            as.data.frame(meta[good]),
            design = reformulate(c(confounders, v)))
        dds <- estimateSizeFactors(dds, type = "poscount")
        dds <- DESeq(dds, parallel = TRUE, quiet = TRUE, fitType = "local")
        res <- lfcShrink(dds, coef = length(resultsNames(dds)),
            res = results(dds))
        res <- as.data.table(res)
        res$genus <- colnames(counts)
        res$variable <- v
        res$n_test <- sum(good)
        tests <- rbind(tests, res)
    }
} else {
    cat("Old test file found. Skipping most tests...\n")
}

# As final variable use diabetes state, this way you can play around with the
# DeSeq object afterwards
good <- as.integer(meta$status) < 7
meta[, status := as.factor(status)]
dds <- DESeqDataSetFromMatrix(t(counts[good, ]), as.data.frame(meta[good]),
                              design = ~ gender + status)
dds <- estimateSizeFactors(dds, type = "poscount")
dds <- DESeq(dds, parallel = TRUE)

combinations <- combn(5:1, 2)
multi <- apply(combinations, 2, function(comb) {
    name <- paste0("status_", comb[1], "_vs_", comb[2])
    res <- results(dds, contrast = c("status", comb[1], comb[2]))
    res <- lfcShrink(dds, contrast = c("status", comb[1], comb[2]), res = res)
    res <- as.data.table(res)
    res$genus <- colnames(counts)
    res$variable <- name
    res$n_test <- min(meta[, table(status)[comb]])
    res
})
multi <- rbindlist(multi)
# Remove genus with a small baseMean to avoid a bimodal pval distribution
multi <- rbind(multi, tests)[baseMean >= min_abundance]
multi[, padj := p.adjust(pvalue, method="fdr")]

if (!skip_all)
    fwrite(multi[order(padj, variable)], "../data/tests_family.csv")
