library(mbtools)

ps <- readRDS("../data/taxonomy.rds")
meta <- readxl::read_excel("../data/clinical_new.xlsx", sheet = "clean")
annotation <- readxl::read_excel("../data/clinical_new.xlsx",
                                 sheet = "annotations_clean")
meta <- meta %>% types(annotation) %>% as.data.frame()
meta <- meta[meta$diabetes_status < 7, ]

# remove duplicates and missing
s <- tstrsplit(sample_names(ps), "_")[[1]]  # sample ids
dupes <- s %in% (which(table(s) != 1) %>% names())
include <- !dupes & (s %in% meta$id)
ps <- prune_samples(include, ps)
sample_names(ps) <- s[include]
meta <- meta[meta$id %in% s[include], ]
rownames(meta) <- meta$id
sample_data(ps) <- meta[s[include], ]
richness <- estimate_richness(ps)
sample_data(ps)$richness <- richness$Chao1
saveRDS(ps, "../data/taxonomy_clean.rds")

# Save the same info in csv files as well
genera <- as.matrix(taxa_count(ps, "genus"))
write.csv(genera, "../data/genera.csv")
write.csv(t(otu_table(ps)), "../data/variants.csv")
write.csv(tax_table(ps), "../data/taxa.csv")

# Save prefiltered and normalized abundance tables

## Filter treated samples
ps <- subset_samples(ps, diabetes_status < 6 & metformin == 0)

genera <- as.matrix(taxa_count(ps, "genus"))
genera <- genera[, colMeans(genera) > 10]
fraction <- apply(genera, 2, function(x) sum(x > 0) / length(x))
genera <- mbtools::normalize(genera[, fraction > 0.1])
write.csv(genera, "../data/counts_norm_genus.csv")

variants <- as(otu_table(ps), "matrix")
variants <- variants[, colMeans(variants) > 10]
fraction <- apply(variants, 2, function(x) sum(x > 0) / length(x))
variants <- mbtools::normalize(variants[, fraction > 0.1])
write.csv(variants, "../data/counts_norm_variant.csv")
