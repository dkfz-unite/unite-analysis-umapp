feature_selection <- function(data, method, n_features) {
  if (is.null(method)) {
    return(data)
  }
  method <- tolower(method)
  if (n_features >= ncol(data)) {
    return(data)
  } else if (method == "variance") {
    feature_variances <- apply(data, 2, var, na.rm = TRUE)
    selected_features <- names(sort(feature_variances, decreasing = TRUE))[1:n_features]
  } else {
    stop(paste("Unsupported feature selection method:", method))
  }
    return (data[, selected_features, drop = FALSE])
}