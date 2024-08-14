## [0.12.0] - 2024-08-14

## Fixed

- fixed a bug in `process_deg_results` where scores, pvals, and pvals_adj were sorted by different means and so
    the records were not in the same order

## [0.11.0] - 2024-08-13

## Fixed

- Added `categorical_covariate_keys` and `continuous_covariate_keys` to `calc_defs()` so that the covariate keys 
    actually get passed on
    - Also added tests to `perform_checks` to make sure or warn about missing or inappropriate keys
- Fixed a bug in `process_anndata` where passing the name of a single layer as `str` instead of a `list[str]` threw an
    error

## [0.10.0] - 2024-08-01

### Changed

- renamed the main function from `scvi_de` to `calc_defs`
- only run `diff_expr_test` in `calc_defs` if `return_df` is True
- handle `save_model_path` being passed to `calc_defs` as a string
- use `muon.pp.filter_var` to subset the anndata object in `create_model` instead of trying 
    to subset the object, avoiding problems with "is it or isn't it a view?"

### Fixed

- use `pd.Series.intersection` in `diff_expr_test` as a function instead of trying to 
    subset it... since it is a function.
- Convert the anndata.X from a `scipy.sparse.csr_matrix` to a `sparse.COO` in `is_integer_array`
    so we can perform math on it, since csr_matrices don't (did they ever?)
    

## [0.9.0] - 2024-07-31

### Added

- pre-commit, nox, mypy, coverage

### Changed

- Removed remnants of `rye`
- Removed `towncrier`

### Fixed

- Use `sparse` to convert `scipy.sparse.csr_matrix` objects to `sparse.COO` before check if they are all integers

## [0.8.0](https://github.com/milescsmith/scvi_de/tree/0.8.0) - 2024-04-29

### Changed

- Now we use `gene-dispersion` and `protein-dispersion` in the correct places (esp different for TOTALVI)

## [0.7.0](https://github.com/milescsmith/scvi_de/tree/0.7.0) - 2024-04-29

### Added

- Ability to control for background proteins (i.e. isotype controls) in CITE-seq data.
    - Note: I'm not entirely sure I'm doing it correctly - copying the counts to obs, then telling TotalVI to use those
    as continuous covariates

### Changed

- DEG results are now filtered to make sure they match what is actually present in the anndata object
- Refactored `process_deg_results`

### Fixed

- Reusing a model no longer tries to reprocess the not-passed anndata/mudata object
- Results are correctly written to `rank_genes_groups` for each modality


### Unchanged

- Documentation remains a mess

## [0.6.0](https://github.com/milescsmith/scvi_de/tree/0.6.0) - 2024-04-26

### Added

- Now works with MuData objects
- A logging submodule, though not using it much at the moment

### Changed

- More properly handle using joint RNA and protein data with TotalVI.
- Changed how differential expression testing works, now using the built-in ability of
    scvi.model.MODEL.differential_expression to test multiple groups. This will probably
    futz with adjusted p_values
- Switch from using rye as as package manager to pdm
- if adata.X is not an integer array, scvi_de now looks to see if adata.raw.X is and copies that
    over if it is
- Remove towncrier

### Refactored

- Split functions up and changed a LOT of when and how things are tested

NOTE: parameter documentation and typing is now a mess

## [0.5.0](https://github.com/milescsmith/scvi_de/tree/0.5.0) - 2024-04-19

### Changed

- Make sure that the dispersion argument to scvi_de is valid

### Fixed

- Did you know that `scipy.sparse.csr_matrix.todense()` and `scipy.sparse.csr_matrix.toarray()`
produce different outputs? Yes? Well apparently I'm an idiot that doesn't notice what he is typing.
Fixed the issue with raw counts being unusable in a layer instead of X because a np.matrix isn't capable of blah blah.


## [0.4.1](https://github.com/milescsmith/scvi_de/tree/0.4.1) - 2024-04-16


### Changed

- Changed the default values for return_model and return_df to True


## [0.4.0](https://github.com/milescsmith/scvi_de/tree/0.4.0) - 2024-04-16


### Added

- Added several points where potential errors are check before expensive calculations are performed

### Fixed

- Performing DEG on protein data actually works now.


## [0.3.4](https://github.com/milescsmith/scvi_de/tree/0.3.4) - 2024-04-16


### Fixed

- Actually pass the type of data from scvi_de to create_model.


## [0.3.3](https://github.com/milescsmith/scvi_de/tree/0.3.3) - 2024-04-16


### Fixed

- Fixed how the layer/X and scVI work (no more tensors full of nan).


## [0.3.2](https://github.com/milescsmith/scvi_de/tree/0.3.2) - 2024-04-15

### Fixed

- Fixed how the layer/X and scVI work (no more tensors full of nan).

## [0.3.1](https://github.com/milescsmith/scvi_de/tree/0.3.1) - 2024-04-15

### Fixed

- Account for if adata.X is a sparse matrix when testing for raw integer counts ([#3](https://github.com/milescsmith/scvi_de/issues/3))

## [0.3.0](https://github.com/milescsmith/scvi_de/tree/0.3.0) - 2024-04-15

### Added

- Added ability to return the model produced and to use a model as input, allowing one to bypass lengthy model creation when changing comparison options
- Allow use of scvi.model.TOTALVI in addition to scvi.model.SCVI ([#20240415](https://github.com/milescsmith/scvi_de/issues/20240415))

## [0.2.0](https://github.com/milescsmith/scvi_de/tree/0.2.0) - 2024-04-15

- add return of scvi model option to scvi_de; add inplace argument ([#20240411](https://github.com/milescsmith/scvi_de/issues/20240411))
