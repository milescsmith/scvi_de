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