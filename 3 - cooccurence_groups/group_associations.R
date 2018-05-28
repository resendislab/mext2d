# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Scripts for running association tests on the genus level

devtools::load_all("../../microbiome/mbtools")
library(magrittr)
library(DESeq2)

ps <- readRDS("../data/taxonomy_clean.rds")
ps <- subset_samples(ps, diabetes_status < 6 & metformin == 0)
variables <- names(sample_data(ps))
exclude <- grepl("_6months", variables) | grepl("_12months", variables) |
           variables %in% c("id", "treatment_group", "metformin")

groups <- fread("sig_healthy_groups.txt")
summarize_groups <- function(taxa, ps) {
    in_group <- taxa_names(ps)[which(tax_table(ps)[, "Genus"] %in% taxa)]
    res <- rowSums(otu_table(ps)[, in_group])
    names(res) <- rownames(otu_table(ps))
    return(res)
}
counts <- groups[, list(counts=list(summarize_groups(taxa, ps))),
                 by = "group_id"]
counts <- do.call(cbind, counts$counts)
gids <- unique(groups$group_id)
counts <- otu_table(cbind(counts), taxa_are_rows = FALSE)
colnames(counts) <- gids
taxa <- matrix(gids)
colnames(taxa) <- "Genus"
rownames(taxa) <- taxa[, 1]
ps <- phyloseq(counts, tax_table(taxa), sample_data(ps))
ps <- subset_samples(ps, apply(otu_table(ps), 1, var) > 0)

tests <- association(ps, variables = variables[!exclude],
                     confounders = c("gender"))

fwrite(tests[order(padj, variable)], "../data/tests_groups.csv")
