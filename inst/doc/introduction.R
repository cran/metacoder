## ----home_setup, echo=FALSE, warning=FALSE, message=FALSE---------------------
options(width = 80)
set.seed(1)
# Knitr
library(knitr)
library(grid)
opts_chunk$set(dev = 'png', fig.width = 7, fig.height = 7, warning = FALSE, message = FALSE)

## ----home_parse, warning=FALSE, message=FALSE---------------------------------
# Load the package
library(metacoder)
# Load the input FASTA file
seqs <- seqinr::read.fasta(system.file("extdata",
                                       "mothur_16S_training_subset.fasta.gz",
                                       package = "metacoder"))
# Print an example of the sequence headers
cat(names(seqs)[1])
# Extract the taxonomic information of the sequences
data <- extract_taxonomy(seqs, regex = "^(.*)\\t(.*)",
                         key = c(id = "obs_info", "class"),
                         class_sep = ";")

## -----------------------------------------------------------------------------
print(data)

## -----------------------------------------------------------------------------
taxon_data(data)

## ----home_plot_1--------------------------------------------------------------
heat_tree(data, node_size = n_obs, node_label = name, node_color = n_obs)

## ----home_plot_3--------------------------------------------------------------
set.seed(1)
heat_tree(filter_taxa(data, name == "Archaea", subtaxa = TRUE),
          node_size = n_obs, node_label = name, 
          node_color = n_obs, layout = "fruchterman-reingold")

## ----home_plot_4--------------------------------------------------------------
subsetted <- filter_taxa(data, n_supertaxa > 0)
set.seed(2)
heat_tree(subsetted, node_size = n_obs, node_label = name,
          node_color = n_obs, tree_label = name)

## ----home_plot_5--------------------------------------------------------------
set.seed(1)
filter_taxa(data, n_supertaxa <= 4) %>%
  heat_tree(node_size = n_obs, node_label = name, node_color = n_obs)

## ----home_plot_sample---------------------------------------------------------
set.seed(1)
sampled <- taxonomic_sample(subsetted, max_counts = c("3" = 10, "5" = 1), min_counts = c("5" = 1))
sampled <- filter_taxa(sampled, n_obs > 0, subtaxa = FALSE) 

## ----home_plot_6--------------------------------------------------------------
set.seed(3)
heat_tree(sampled, 
          node_size = n_obs,
          node_label_size = n_obs * ifelse(n_supertaxa == 3, 10, 1),
          edge_size = n_obs, 
          node_label = n_obs,
          node_color = n_obs,
          tree_label = name)

## ----home_plot_7--------------------------------------------------------------
set.seed(6)
sample_n_obs(subsetted, size = 400, taxon_weight = 1 / n_obs, unobserved = FALSE) %>%
  heat_tree(node_size = n_obs, node_label = n_obs, overlap_avoidance = 0.5,
            node_color = n_obs, tree_label = name)

