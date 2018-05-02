## ----home_setup, echo=FALSE, warning=FALSE, message=FALSE-------------------------------
options(width = 90)
set.seed(1)
# Knitr
library(knitr)
library(grid)
opts_chunk$set(dev = 'png', fig.width = 7, fig.height = 7, warning = TRUE,
               message = TRUE)

## ---------------------------------------------------------------------------------------
library(metacoder)
print(hmp_otus)
print(hmp_samples)

## ---------------------------------------------------------------------------------------
obj <- parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                      class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                      class_regex = "^(.+)__(.+)$")

## ---------------------------------------------------------------------------------------
print(obj)

## ---------------------------------------------------------------------------------------
obj$data$tax_data <- zero_low_counts(obj, "tax_data", min_count = 5)

## ---------------------------------------------------------------------------------------
no_reads <- rowSums(obj$data$tax_data[, hmp_samples$sample_id]) == 0
sum(no_reads)

## ---------------------------------------------------------------------------------------
obj <- filter_obs(obj, "tax_data", ! no_reads, drop_taxa = TRUE)
print(obj)

## ---------------------------------------------------------------------------------------
obj$data$tax_data <- calc_obs_props(obj, "tax_data")
print(obj)

## ---------------------------------------------------------------------------------------
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = hmp_samples$sample_id)
print(obj)

## ---------------------------------------------------------------------------------------
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = hmp_samples$body_site)
print(obj)

## ---------------------------------------------------------------------------------------
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Nose, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads")

## ---- eval = FALSE----------------------------------------------------------------------
#  heat_tree(obj,
#            node_label = obj$taxon_names(),
#            node_size = obj$n_obs(),
#            node_color = obj$data$tax_occ$Nose,
#            node_size_axis_label = "OTU count",
#            node_color_axis_label = "Samples with reads")

## ---- warning = FALSE-------------------------------------------------------------------
obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund",
                                      cols = hmp_samples$sample_id,
                                      groups = hmp_samples$sex)
print(obj$data$diff_table)

## ---------------------------------------------------------------------------------------
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          edge_color_interval = c(-2, 2),
          node_color_range = c("cyan", "gray", "tan"),
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median proportions")

## ---------------------------------------------------------------------------------------
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
                                               method = "fdr")

## ---------------------------------------------------------------------------------------
hist(obj$data$diff_table$wilcox_p_value) 

## ---- warning = FALSE-------------------------------------------------------------------
obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund",
                                      cols = hmp_samples$sample_id,
                                      groups = hmp_samples$body_site)
print(obj$data$diff_table)

## ---------------------------------------------------------------------------------------
heat_tree_matrix(obj,
                 dataset = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions")

