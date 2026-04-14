library(limma)
library(readr)
library(dplyr)
library(tibble)
library(jsonlite)
library(umap)
source(file.path(getwd(), "helpers", "preprocessing.r"))
source(file.path(getwd(), "helpers", "feature_selection.r"))
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
metadata_file <- file.path(workdir, "metadata.tsv")
results_file <- file.path(workdir, "results.tsv")
data_file <- file.path(workdir, "data.tsv")
options_file <- file.path(workdir, "options.json")

# Read data and metadata
data <- read_tsv(data_file)
metadata <- read_tsv(metadata_file)
options <- fromJSON(options_file)

### preprocess options
# set umap random_state to NA_integer_ if not provided or empty, otherwise convert to integer
rs <- get_required(options, "umap_random_state")
if (is.null(rs) || length(rs) == 0) {
rs <- NA_integer_
} else {
rs <- as.integer(rs)
}
options <- replace_required(options, "umap_random_state", rs)

# if required_min_fraction_one_class is TRUE, check that condition column contains no missing values
if (get_required(options, "require_min_fraction_one_class")) {
  if (anyNA(metadata$condition)) {
    stop("require_min_fraction_one_class is TRUE but condition column contains missing values")
  }
}



print(names(options))
# Preprocess data (log, normalize, impute, batch correct)
data_matrix <- as.data.frame(data[,-1]) # remove feature column for processing
rownames(data_matrix) <- data[[1]] # set feature names as rownames

metadata_matrix <- as.data.frame(metadata[, -1]) # assuming first column is sample names
rownames(metadata_matrix) <- metadata[[1]] # set sample names as rownames

# reorder metadata rows to match data matrix columns
metadata_matrix <- metadata_matrix[match(colnames(data_matrix), rownames(metadata_matrix)), ,drop = FALSE]

# transpose data_matrix as proteomic_data_preprocessing expects samples as rows, features as columns
data_matrix <- t(data_matrix)

# preprocess data
processed_data <- preprocess_data(data=data_matrix, 
                              batch_vector=get_required(metadata_matrix, "batch"),
                              class_labels=get_required(metadata_matrix, "condition"), 
                              options = options)


# do feature selection
processed_data <- feature_selection(processed_data,
                method = get_required(options, "feature_selection_method"),
                n_features = get_required(options, "feature_selection_n_features"))


### run initial dimensionality reduction with pca to reduce computational burden
pca_result <- prcomp(processed_data, center = TRUE, scale. = TRUE)
n_comps_to_retain <- min(get_required(options, "umap_n_principal_components"), ncol(pca_result$x))
pca_scores <- pca_result$x[, 1:n_comps_to_retain, drop = FALSE]


# unpack umap configuration options
umap_config <- umap.defaults
umap_config <- replace_required(umap_config, "n_neighbors", min(get_required(options, "umap_n_neighbors"), nrow(pca_scores) - 1))
umap_config <- replace_required(umap_config, "metric", tolower(get_required(options, "umap_metric")))
umap_config <- replace_required(umap_config, "random_state", get_required(options, "umap_random_state"))
umap_config <- replace_required(umap_config, "min_dist", get_required(options, "umap_min_dist"))
print(umap_config)
# run umap on pca scores
umap_result <- umap(pca_scores, config = umap_config)

# write umap coordinates to file
umap_coords <- as.data.frame(umap_result$layout)
umap_coords$sample <- rownames(pca_scores)
umap_coords <- dplyr::select(umap_coords, sample, dplyr::everything())
write_tsv(umap_coords, file.path(workdir, "results.tsv"))

