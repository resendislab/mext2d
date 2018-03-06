library(missForest)
library(phyloseq)

ps <- readRDS("../data/taxonomy_clean.rds")
meta <- as(sample_data(ps), "data.frame")

variables <- names(sample_data(ps))
exclude <- grepl("months", variables) | variables == "id"
meta <- meta[, !exclude]
meta <- missForest(meta)$ximp
