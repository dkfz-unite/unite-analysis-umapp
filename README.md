### Data format
Data file `{proc}/data.tsv` should be in the following format:
```tsv
feature sample1 sample2 sample3 sample4
protein1 10 20 30 40
protein2 15 25 35 45
protein3 20 30 40 50
protein4 25 35 45 55
```  

Where:
- `feature` - identifier of the feature (protein). Should be first column.
- `sample1`, `sample2`, `sample3`, `sample4` - names of the samples.
- Values in the table are raw protein intensity values (not normalissed, not filtered).

### Metadata format
Metadata file `{proc}/metadata.tsv` should be in the following format:
```tsv
sample batch metadata1 metadata2
sample1 1 8 k
sample2 2 4 k
sample3 1 2 k
sample4 2 1 k
```

Where:
- `sample` - name of the sample. Should be first column.
- `batch` - optional batch variable for batch correction. If there is no batch variable, this column should be present but empty
- `metadata1` and `metadata2` represent optional columns that can be present, but will have no effect on the analysis.

### Options configuration
The preprocessing of the data is configurable.

- `normalization_method`: ["median", "quantile"] data are first log2 normalized then either median centered (if "median") or quantile normalized using `preprocessCore`.
- `normalization_log_offset`:[float>0] log2 transformation applies an offset `log2(x + offset)` 
- `imputation method`: ["mindet", "minprob"] missing values (NA or 0) are imputed using either:
  - "mindet" algorithm where each missing value is imputed with the corresponding "minimum" (actually 1st percentile) of intensity values for its sample.
  - "minprob" algorithm where imputed values are sampled from a Gaussian distribution, centered on the "minimum" (see above). The standard deviation of this distribution is set to the median of the standard deviations of the distributions of all proteins.
- `stratify_imputation_by_batch`: [true, false] if true the imputation will be done separately for each batch. In the event that any proteins have less than three non-missing samples in any batch, it will fall back to non-stratified processing.
- `batch_correction_method`: ["combat", "limma", null] if there is a batch variable will perform batch correction with `comBat` (par.prior=True) function of the `sva` package or `removeBatchEffects` from `limma`. 
- `min_non_missing_fraction`:[float (0<x<=1.)] the minimum proportion of non-missing values of a protein required for it be retained for analysis
- `require_min_fraction_one_class`: [true,false] if true a protein must exceed `min_non_missing_fraction` in only one class, to retained.
- `umap_n_neighbors` [int 2<=x<n_samples] the number opf neighbors umap will use to build the graph
- `umap_metric` [string literal ["euclidean", "manhattan"]] the distance metric to use for determining similarity among the data points.
- `umap_random_state` [int >=0 or null] seed to ensure reproducibility.
- `umap_min_dist` [float > 0.] determines how close together points can be in the final embedding.
- `umap_n_principal_components` [int >= 1] for high dimensional data the UMAP algorithm becomes very expensive. This parameter specifies the number of principal components on which to base the final embedding.
- `feature_selection_method` [string literal ["variance", "none"]] optionally perform feature selection prior to PCA. If "variance" this will select the top n features with the highest variance.
- `feature_selection_n_features` [int >= 1] the number of features to select if feature_selection_method is not "none".
```json
{
  "normalization_method": "median",
  "normalization_log_offset": 1,
  "imputation_method": "mindet",
  "stratify_imputation_by_batch": false,
  "batch_correction_method": null,
  "min_non_missing_fraction": 0.5,
  "require_min_fraction_one_class": false,

  "feature_selection_method": "variance",
  "feature_selection_n_features": 1000,

  "umap_n_neighbors": 15,
  "umap_metric": "euclidean",
  "umap_random_state": null,
  "umap_min_dist": 0.1,
  "umap_n_principal_components": 50
}
```
