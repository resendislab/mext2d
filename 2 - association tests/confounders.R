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

vars <- c("diabetes_status", "glucose_0", "glucose_120", "auc_glucose",
          "glycated_haemoglobin", "insulin_0", "auc_insulin")
confounders <- c("gender", "bmi", "waist_hip_ratio", "percent_body_fat",
                 "visceral_fat")
diabetes_no_obesity <- association(ps, variables = vars,
                                   confounders = confounders)

vars <- c("gender", "bmi", "waist_hip_ratio", "percent_body_fat",
          "visceral_fat", "bmi_status")
confounders <- c("gender", "diabetes_status", "glucose_0", "glucose_120",
                 "auc_glucose", "glycated_haemoglobin", "insulin_0",
                 "auc_insulin")
obesity_no_diabetes <- association(ps, variables = vars,
                                   confounders = confounders)

vars <- c("systolic_pressure", "diastolic_pressure", "mean_blood_pressure",
          "pulse_pressure", "cholesterol")
confounders <- c("gender", "diabetes_status", "glucose_0", "glucose_120",
                 "auc_glucose", "glycated_haemoglobin", "insulin_0",
                 "auc_insulin")
cardio_no_diabetes <- association(ps, variables = vars,
                                  confounders = confounders)
