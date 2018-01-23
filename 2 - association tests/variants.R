# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Scripts for running association tests on the genus level

library(mbtools)
library(magrittr)
library(DESeq2)

zero_fraction <- 0.1
min_abundance <- 8
skip_all <- file.exists("../data/tests_species.csv")

ps <- readRDS("../data/taxonomy.rds")
taxa <- taxa_count(ps, lev=NA)
setkey(taxa, "taxa")
counts <- as.matrix(taxa)

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
fraction_exclude <- colSums(counts >= 1) / nrow(counts) < zero_fraction
cat("removed variants due to missing reads:",
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
        res$variant <- colnames(counts)
        res$variable <- v
        res$n_test <- sum(good)
        tests <- rbind(tests, res)
    }
} else {
    cat("Old test file found. Skipping most tests...\n")
}

# As final variable use diabetes state, this way you can play around with the
# DeSeq object afterwards
good <- as.integer(meta$status) < 6
meta[, status := as.factor(status)]
dds <- DESeqDataSetFromMatrix(t(counts[good, ]), as.data.frame(meta[good]),
                              design = ~ gender + status)
dds <- estimateSizeFactors(dds, type = "poscount")
dds <- DESeq(dds, parallel = TRUE)

multi <- NULL
for (level in 2:5) {
    name <- paste0("status_", level, "_vs_1")
    res <- results(dds, name = name)
    res <- lfcShrink(dds, coef = which(resultsNames(dds) == name), res = res)
    res <- as.data.table(res)
    res$variant <- colnames(counts)
    res$variable <- name
    res$n_test <- nrow(meta[status == level])
    multi <- rbind(multi, res)
}
# Remove genus with a small baseMean to avoid a bimodal pval distribution
multi <- rbind(multi, tests)[baseMean >= min_abundance]
multi[, padj := p.adjust(pvalue)]
variant_map <- taxa[, .(species = unique(species)), by="taxa"]
setkey(variant_map, "taxa")
multi[, species := variant_map[variant, species]]

if (!skip_all)
    fwrite(multi[order(padj, variable)], "../data/tests_variants.csv")
