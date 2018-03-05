devtools::load_all("../../microbiome/mbtools")

ps <- readRDS("../data/taxonomy.rds")
meta <- readxl::read_excel("../data/clinical_new.xlsx", sheet = "clean")
annotation <- readxl::read_excel("../data/clinical_new.xlsx", sheet = "annotations_clean")
meta <- meta %>% types(annotation) %>% standardize() %>% as.data.frame()
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
saveRDS(ps, "../data/taxonomy_clean.rds")

