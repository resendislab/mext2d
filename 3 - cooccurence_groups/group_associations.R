# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Scripts for running association tests on the genus level

devtools::load_all("../../microbiome/mbtools")
library(magrittr)
library(DESeq2)

ps <- readRDS("../data/taxonomy_clean.rds")
ps <- subset_samples(ps, diabetes_status < 6 && metformin == 0)
variables <- names(sample_data(ps))
exclude <- grepl("_6months", variables) | grepl("_12months", variables) |
           variables == "id"

groups <- fread("sig_cooccurence_groups.txt")
groups <- groups[, .(correlated=list(correlated_taxa)),
                 by = c("taxa", "correlation", "selection_method")]
groups$name <- paste0("group_", 1:nrow(groups))
summarize_groups <- function(taxa, ps) {
    taxa <- unlist(taxa)
    in_group <- taxa_names(ps)[which(tax_table(ps)[, "Genus"] %in% taxa)]
    res <- rowSums(otu_table(ps)[, in_group])
    names(res) <- rownames(otu_table(ps))
    return(res)
}
counts <- groups[, list(counts=list(summarize_groups(correlated, ps))),
                 by = "name"]
counts <- do.call(cbind, counts$counts)
colnames(counts) <- groups$name
counts <- otu_table(cbind(counts), taxa_are_rows = FALSE)
taxa <- matrix(groups$name)
colnames(taxa) <- "Genus"
rownames(taxa) <- taxa[, 1]
ps <- phyloseq(counts, tax_table(taxa), sample_data(ps))

tests <- association(ps, variables = variables[!exclude],
                     confounders = c("gender"))
tests <- tests[groups, on = c(genus="name")]
tests[, genus := NULL]

fwrite(tests[order(padj, variable)], "../data/tests_groups.csv")
