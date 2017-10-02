#!/usr/bin/env Rscript

# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0.

library(mbtools)
library(stringr)
library(data.table)


annotate_files <- function(path) {
    fwd <- list.files(path, pattern="R1.+\\.fastq.gz", full.names = TRUE,
                      recursive = TRUE)
    ids_fwd <- str_match(fwd, "(Placa.+)/(.+)_S")
    bwd <- list.files(path, pattern = "R2.+\\.fastq.gz", full.names = TRUE,
                      recursive = TRUE)
    ids_rev <- str_match(bwd, "(Placa.+)/(.+)_S")
    files_fwd <- data.table(run = ids_fwd[, 2],
                            sample = ids_fwd[, 3], forward = fwd)
    files_rev <- data.table(run = ids_rev[, 2],
                            sample = ids_rev[, 3], reverse = bwd)
    files <- files_fwd[files_rev, on = .(plate, run, sample)]

    return(files)
}

dada_errors <- function(samples) {
    fwd_err <- learnErrors(samples$forward, nreads=2e6, multithread=TRUE,
                           randomize=TRUE)
    rev_err <- learnErrors(samples$reverse, nreads=2e6, multithread=TRUE,
                           randomize=TRUE)

    return(list(forward=list(fwd_err), reverse=list(rev_err)))
}

samples <- fread("samples.csv")
samples <- add_files("filtered", samples)

if (!file.exists("errors.rds")) {
    err <- samples[, dada_errors(.SD), by="run"]
    saveRDS(err, "errors.rds")
} else {
    err <- readRDS("errors.rds")
}
setkey(err, run)

for(r in unique(samples$run)) {
    cat("Processing run", r, "...\n")
    filename <- paste0("merged_", r, ".rds")
    if (file.exists(filename)) next
    s <- samples[run == r]
    derep_fwd <- derepFastq(s$forward)
    derep_rev <- derepFastq(s$reverse)
    dd_fwd <- dada(derep_fwd, err = err[r, forward][[1]], multithread = TRUE)
    dd_rev <- dada(derep_rev, err = err[r, reverse][[1]], multithread = TRUE)
    merged <- mergePairs(dd_fwd, derep_fwd, dd_rev, derep_rev)
    saveRDS(merged, filename)
}

merged <- lapply(list.files(pattern = "merged_"),
                 function(path) makeSequenceTable(readRDS(path)))
seqtab <- do.call(mergeSequenceTables, merged)
seqtab <- removeBimeraDenovo(seqtab, multithread = TRUE)
saveRDS(seqtab, "seqtab.rds")

taxa <- assignTaxonomy(seqtab, "../silva_nr_v128_train_set.fa.gz",
                       multithread = TRUE)
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               tax_table(taxa))
saveRDS(ps, "taxonomy.rds")
