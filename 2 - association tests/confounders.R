# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# Some association tests with additional confounders

devtools::load_all("../../microbiome/mbtools")
library(magrittr)
library(DESeq2)

ps <- readRDS("../data/taxonomy_clean.rds")
ps <- subset_samples(ps, diabetes_status < 6 & metformin == 0)
variables <- names(sample_data(ps))
exclude <- grepl("_6months", variables) | grepl("_12months", variables) |
           variables %in% c("id", "treatment_group", "metformin")

diabetes <- c("diabetes_status", "auc_glucose", "glycated_haemoglobin",
              "auc_insulin", "num_risk_factors")
obesity <- c("bmi", "waist_hip_ratio", "visceral_fat", "percent_body_fat")
cardio <- c("systolic_pressure", "diastolic_pressure", "cholesterol")

tests <- list()
tests$diabetes <- association(ps, variables = diabetes,
                              confounders = c("gender", obesity, cardio))
tests$diabetes[, "response" := "diabetes"]

tests$obesity <- association(ps, variables = obesity,
                             confounders = c("gender", diabetes, cardio))
tests$obesity[, "response" := "obesity"]

tests$cardio <- association(ps, variables = cardio,
                            confounders = c("gender", diabetes, obesity))
tests$cardio[, "response" := "cardio"]

tests <- rbindlist(tests)
fwrite(tests, "../data/tests_confounders.csv")
