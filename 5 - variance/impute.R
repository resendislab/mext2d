library(missForest)
library(phyloseq)

ps <- readRDS("../data/taxonomy_clean.rds")
meta <- as(sample_data(ps), "data.frame")

variables <- names(sample_data(ps))
exclude <- grepl("months", variables) | variables %in% c("id", "treatment_group")
meta <- meta[, !exclude]
meta <- missForest(meta)$ximp
