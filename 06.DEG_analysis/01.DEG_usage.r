library(dplyr)
library(tidyverse)

# load data
base_dir <- "/disk2/cai113/data/stateTrans"
count_file <- file.path(base_dir, "04.standardization", "01.counts_raw.txt")
count_mat <- as.matrix(read.delim(count_file, row.names = 1, check.names = FALSE))
load(file.path(base_dir, "05.model", "model_results.rdata"))
# load scripts
source(file.path(base_dir, "06.DEG_analysis", "00.DEG_DESeq2_script.r"))

# create sample info
sample <- c(model_results$x_wt, model_results$x_ko) %>% as.data.frame()
sample$name <- rownames(sample)
colnames(sample)[1] <- "PC"
# cluster by kmeans results
kms <- model_results$critical_points$kmeans_result$cluster
sample$kms <- ifelse(is.na(kms[sample$name]), 0, kms[sample$name])
# cluster by critical points. however, due to lack of samples, and the similar pattern with kmeans, I decide to use kmeans to refer critical points.
sample$cp <- ifelse(sample$kms == 0, "c1_star", ifelse(sample$kms == 3, "c1", ifelse(sample$kms == 2,"c2", "c3")))

# run DEG analysis batch
contrast_list <- list(c("cp", "c1", "c1_star"),
                      c("cp", "c2", "c1_star"),
                      c("cp", "c3", "c1_star"),
                      c("cp", "c2", "c1"),
                      c("cp", "c3", "c1"),
                      c("cp", "c3", "c2"))
batch_results <- batch_deseq2_analysis(
    count_matrix = count_mat,
    sample_info = sample,
    group_col = "cp",
    contrast_list = contrast_list,
    output_dir = file.path(base_dir, "06.DEG_pathway", "DEG_results"),
    create_subdir = TRUE)