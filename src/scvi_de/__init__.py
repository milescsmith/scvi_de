import warnings
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Literal

import muon as mu
import numpy as np
import numpy.typing as npt
import pandas as pd
import scvi
import sparse
from anndata import AnnData
from loguru import logger
from mudata import MuData
from scipy.sparse import csr_matrix, issparse

from scvi_de.logging import init_logger

init_logger(verbose=3)

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"


# TODO: add logging
# TODO: pass more options on
# TODO: add an option to save the model IMMEDIATELY after creation
# TODO: FUCKING HELL WRITE SOME UNIT TESTS!
def calc_defs(
    adata: AnnData | MuData | None = None,
    model: scvi.model.SCVI | scvi.model.TOTALVI | None = None,
    save_model_path: Path | str | None = None,
    overwrite_previous_model: bool = False,
    protein_obsm_key: str = "prot",
    protein_names_uns_key: str = "prot_names",
    groupby: str = "leiden",
    layers: str | dict[str, str] | None = None,
    protein_layer: str | None = None,
    modality: Literal["rna", "prot"] = "rna",
    mudata_protein_modality: str = "prot",
    mudata_rna_modality: str = "rna",
    compare_group: str | None = None,
    reference_group: str | None = None,
    batch_correct: bool = False,
    batch_key: str | None = None,
    categorical_covariate_keys: str | list[str] | None = None,
    continuous_covariate_keys: str | list[str] | None = None,
    gene_dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
    protein_dispersion: Literal["protein", "protein-batch", "protein-label"] = "protein",
    gene_likelihood: Literal["nb", "zinb", "poisson"] = "nb",
    remove_outliers: bool = False,
    return_df: bool = True,
    return_model: bool = True,
    lfc_use: Literal["lfc_mean", "lfc_median", "lfc_max", "lfc_min"] = "lfc_mean",
    use_isotype_controls: bool = True,
    isotype_pattern: str = "Isotype",
    num_devices: str | int = "auto",
    **kwargs,
) -> dict[str, pd.DataFrame | AnnData | scvi.model.SCVI | scvi.model.TOTALVI] | None:
    """Perform differential expression analysis using the scVI suite of tools

    Parameters
    -------
    adata : :class:`~anndata.AnnData` | :class:`~MuData`
        Object to model and perform differential expression on. Either an AnnData/MuData or model object need to be passed;
        if both are passed, AnnData will be ignored.
    model : :class:`~scvi.model.SCVI`, :class:`~scvi.model.TOTALVI`
        A previously trained model.  If both an AnnData and model object are passed, the AnnData object will be ignored
        and the model will be used.
    save_model_path : :class:`~pathlib.Path`
        Path to a folder in which to save the model and associated anndata object. This value will be used to create a
        folder within which will be placed 'model.pt' and 'adata.h5ad' files. This WILL overwrite any e
    protein_obsm_key: str | None, default=None

    protein_names_uns_key: str | None, default=None

    groupby : str, default="leiden"
        Column in `obs` to use when performing comparisons.
    layer : str
        Optional. Key for `adata.layers` where raw count data is stored. If not provided, the counts in `adata.X` will
        be used
    modality : Literal["rna", "prot"], default="rna"
        Type of data to be processed.  Currently, only "rna" (for RNA-seq) and "prot" (for CITE-seq) are allowed.
    compare_group : str
        Category in `adata.obs[groupby]` to use as the query for comparison of two groups. Positive fold
        changes indicate higher expression in this group.
    reference_group : str
        Category in `adata.obs[groupby]` to use as a reference for comparison of two groups. Negative
        fold changes indicate higher expression in this group.
    batch_correct : bool, default=False
        Should batch effects be corrected for in DE inference?
    batch_key : str, default="batch"
        Column in `adata.obs` to use to indicate batches.
    dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell", "protein", "protein-batch", "protein-label"]
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
        * ``'protein'`` - dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - dispersion can differ between different labels NOT TESTED

        If not passed, either 'gene' or 'protein' will be selected based on the type of data.
    gene_likelihood : Literal["nb", "zinb", "poisson"], default="nb"
        One of the following:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution

    remove_outliers : bool, default=False
        Whether to filter outlier cells with `:meth:~scvi.model.base.DifferentialComputation.filter_outlier_cells()`
    return_df : bool, default=True
        Should a dataframe of results be returned?
    return_model : bool, default=True
        Should the model be returned?
    lfc_use : Literal["lfc_mean", "lfc_median", "lfc_max", "lfc_min"], default="lfc_mean"
        Which log fold change to use.

    Return
    ------
        dict of {str: pd.DataFrame | scvi.model.SCVI | scvi.model.TOTALVI}
        Original object is modified in place, with results placed in adata.uns["rank_genes_groups"]
    """

    # several check to make sure we don't waste a lot of time on something that will later fail
    logger.debug("performing checks")
    if isinstance(save_model_path, str):
        save_model_path = Path(save_model_path)

    perform_checks(
        adata,
        model,
        save_model_path,
        overwrite_previous_model,
        protein_obsm_key,
        modality,
        mudata_protein_modality,
        mudata_rna_modality,
        gene_dispersion,
        protein_dispersion,
        categorical_covariate_keys,
        continuous_covariate_keys,
    )

    # scvi sometimes takes MuData objects, sometimes not.  I don't think TotalVI does, so we need to reform the data
    # both scVI and TotalVI models assume RNA data in X or the indicated layer so need this in either case
    if not model:
        tmp = (
            process_mudata(
                adata,
                layers,
            )
            if isinstance(adata, MuData)
            else process_anndata(adata, layers)
        )
    else:
        tmp = model.adata

    if adata and not model:
        model = create_model(
            adata=tmp,
            layers=layers,
            batch_key=batch_key,
            gene_dispersion=gene_dispersion,
            protein_dispersion=protein_dispersion,
            gene_likelihood=gene_likelihood,
            modality=modality,
            # num_devices=num_devices,
            use_isotype_controls=use_isotype_controls,
            isotype_pattern=isotype_pattern,
            categorical_covariate_keys=categorical_covariate_keys,
            continuous_covariate_keys=continuous_covariate_keys,
            **kwargs,
        )
    elif (adata is None) and (model is None):
        msg = "Either an AnnData/MuData object or a scVI/TotalVI model are needed and both are missing."
        raise ValueError(msg)

    if save_model_path:
        model.save(save_model_path, overwrite=overwrite_previous_model, save_anndata=True)

    return_dict = {}

    if return_df:
        de_change = diff_expr_test(
            adata=model.adata,
            model=model,
            groupby=groupby,
            current_group=compare_group,
            reference_group=reference_group,
            batch_correct=batch_correct,
            remove_outliers=remove_outliers,
        )

        if isinstance(model.adata, AnnData):
            model.adata.uns["rank_genes_groups"] = process_deg_results(
                df=de_change,
                groupby=groupby,
                layer=layers,
                lfc_use=lfc_use,
            )
        elif isinstance(model.adata, MuData):
            for _ in model.adata.mod:
                if isinstance(de_change, dict):
                    model.adata[_].uns["rank_genes_groups"] = process_deg_results(
                        df=de_change[_], groupby=groupby, layer=layers[_] if layers else None, lfc_use=lfc_use
                    )

        return_dict["deg_df"] = de_change

    model.adata.obsm["X_scVI"] = model.get_latent_representation()

    if return_model:
        return_dict["model"] = model

    return return_dict


def perform_checks(
    adata: AnnData | MuData,
    model: scvi.model.SCVI | scvi.model.TOTALVI,
    save_model_path: Path | None = None,
    overwrite_previous_model: bool = False,
    protein_obsm_key: str | None = None,
    modality: Literal["rna", "prot"] = "rna",
    mudata_protein_modality: str = "prot",
    mudata_rna_modality: str = "rna",
    gene_dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
    protein_dispersion: Literal["protein", "protein-batch", "protein-label"] = "protein",
    categorical_covariate_keys: str | list[str] | None = None,
    continuous_covariate_keys: str | list[str] | None = None,
):
    if save_model_path and (save_model_path.exists() and not overwrite_previous_model):
        msg = (
            f"It appears that a folder exists at {save_model_path}, but `overwrite_previous_model` is `False`. "
            f"Make sure that is correct before long and expensive calculations are carried out!"
        )
        warnings.warn(msg, UserWarning, stacklevel=2)

    if isinstance(adata, MuData):
        for _ in (mudata_protein_modality, mudata_rna_modality):
            if _ not in adata.mod:
                msg = (
                    f"{_} does not appear to be the name of a modality in the passed object. "
                    f"The available modalities are {','.join(adata.mod)}."
                )
                raise ValueError(msg)
    if (protein_obsm_key is None) and ("prot" in modality.lower()) and isinstance(adata, AnnData):
        msg = (
            "A name for the key where protein data in adata.obsm was not provided, but the modality selected is protein. "
            "If assessing protein data, it is recommended to to either use a MuData object or place transcript counts "
            "in adata.X or a layer in adata and the protein counts in adata.obsm[protein_obsm_key]. Results WILL differ."
        )
        warnings.warn(msg, UserWarning, stacklevel=2)
    elif isinstance(adata, AnnData) and ("prot" in modality.lower()):
        msg = (
            "If assessing protein data, it is recommended to to either use a MuData object or place transcript counts "
            "in adata.X or a layer in adata and the protein counts in adata.obsm[protein_obsm_key]. Results WILL differ."
        )
        warnings.warn(msg, UserWarning, stacklevel=2)

    if gene_dispersion not in (
        "gene",
        "gene-batch",
        "gene-label",
        "gene-cell",
    ):
        msg = f"{gene_dispersion} is not a valid option for gene_dispersion"
        raise ValueError(msg)

    if protein_dispersion not in (
        "protein",
        "protein-batch",
        "protein-label",
    ):
        msg = f"{protein_dispersion} is not a valid option for protein_dispersion"
        raise ValueError(msg)
    if adata and model:
        msg = (
            "When both an Anndata/MuData object and scVI model are passed to scvi_de, the existing model is used for "
            "comparisons and is not regenerated. Instead, the differential expression results are written to "
            "the provided Anndata/MuData object in adata.uns['rank_genes_groups']."
        )
        warnings.warn(msg, UserWarning, stacklevel=2)
    if model:
        msg = (
            "Reusing the provided model. If you want to add the results from the current comparisons to your "
            "AnnData/Mudata object, run `adata.uns['rank_genes_groups'] = process_deg_results(degs['deg_df'])` "
            "or the MuData equivalent."
        )
        warnings.warn(msg, UserWarning, stacklevel=2)
    if categorical_covariate_keys and adata:
        if not_found := [_ for _ in categorical_covariate_keys if _ not in adata.obs.columns]:
            msg = f"{not_found.join(', ')} were passed as covariate keys, but they are not present in the obs columns."
            raise ValueError(msg)
        elif non_category := [
            _ for _ in categorical_covariate_keys if not pd.api.types.is_categorical_dtype(adata.obs[_])
        ]:
            msg = (
                f"{non_category.join(', ')} were passed as categorical covariate keys, but their data type is not "
                "categorical"
            )
            warnings.warn(msg, UserWarning, stacklevel=2)
    if continuous_covariate_keys and adata:
        if not_found := [_ for _ in continuous_covariate_keys if _ not in adata.obs.columns]:
            msg = f"{not_found.join(', ')} were passed as covariate keys, but they are not present in the obs columns."
            raise ValueError(msg)
        elif non_numeric := [_ for _ in categorical_covariate_keys if not pd.api.types.is_numeric_dtype(adata.obs[_])]:
            msg = (
                f"{non_numeric.join(', ')} were passed as categorical covariate keys, but their data type is not "
                "categorical"
            )
            warnings.warn(msg, UserWarning, stacklevel=2)


def process_deg_results(
    df: pd.DataFrame,
    groupby: str,
    layer: str | None = None,
    lfc_use="lfc_mean",
):
    """
    Reformat the dataframe(s) returned by `scvi.model.MODEL.differential_expression()` to a form matching
    that returned by `scanpy.tl.rank_genes_groups()`.

    Parameters
    ----------
    de_change : :class:`pd.DataFrame` | dict[str, :class:`pd.DataFrame`]
        A dataframe for a single comparison or dictionary of dataframes for multiple comparisons, as returned by
        `scvi.model.MODEL.differential_expression()`
    lfc_use : :class:`Literal["lfc_mean", "lfc_median", "lfc_std", "lfc_min", "lfc_max"]`, default="lfc_mean"
        Which log-fold change column in the dataframes to use for fold change values

    Return
    ------
    Same as that returned by `scanpy.tl.rank_genes_groups()`
    """

    df["score"] = df.apply(lambda x: x["proba_de"] * x["lfc_mean"] * x["non_zeros_proportion1"], axis=1)
    score_order = {
        i: df.loc[lambda x: x["group1"] == i].sort_values("score", ascending=False).index  # noqa: B023 I don't know why ruff flags this line
        for i in sorted(df["group1"].unique())
    }

    names_rec = pd.DataFrame(
        {i: df[df["group1"] == i].loc[score_order[i], "score"].index.to_list() for i in sorted(df["group1"].unique())}
    ).to_records(index=False)

    scores_rec = pd.DataFrame(
        {i: df[df["group1"] == i].loc[score_order[i], "score"].to_list() for i in sorted(df["group1"].unique())}
    ).to_records(index=False)

    pvals_rec = pd.DataFrame(
        {i: df[df["group1"] == i].loc[score_order[i], "proba_not_de"].to_list() for i in sorted(df["group1"].unique())}
    ).to_records(index=False)

    logfoldchanges_rec = pd.DataFrame(
        {i: df[df["group1"] == i].loc[score_order[i], lfc_use].to_list() for i in sorted(df["group1"].unique())}
    ).to_records(index=False)

    return {
        "params": {
            "groupby": groupby,
            "reference": df["group2"].unique()[0],
            "method": "scvi_de",
            "use_raw": True,
            "layer": layer,
            "lfc_used": lfc_use,
            "corr_method": "benjamini-hochberg",
        },
        "names": names_rec,
        "scores": scores_rec,
        "pvals": pvals_rec,
        "pvals_adj": pvals_rec,
        "logfoldchanges": logfoldchanges_rec,
    }


def create_model(
    adata: AnnData | MuData,
    layers: str | dict[str, str] | None = None,
    modality: Literal["rna", "prot"] = "rna",  # add more options and we validate input arguments
    mudata_rna_modality: str = "rna",
    mudata_protein_modality: str = "prot",
    use_isotype_controls: bool = True,
    isotype_pattern: str = "Isotype",
    batch_key: str | None = None,
    protein_obsm_key: str | None = None,
    protein_names_uns_key: str | None = None,
    gene_dispersion: Literal[
        "gene",
        "gene-batch",
        "gene-label",
        "gene-cell",
    ] = "gene",
    protein_dispersion: Literal[
        "protein",
        "protein-batch",
        "protein-label",
    ] = "protein",
    gene_likelihood: Literal["nb", "zinb", "poisson"] = "nb",
    size_factor_key: str | list[str] | None = None,
    categorical_covariate_keys: str | list[str] | None = None,
    continuous_covariate_keys: str | list[str] | None = None,
    save_layer: str = "counts",
    check_val_every_n_epoch: int = 1,
    max_epochs: int = 400,
    early_stopping: bool = True,
    early_stopping_patience: int = 20,
    early_stopping_monitor: Literal[
        "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
    ] = "elbo_validation",
    num_devices: str | int = "auto",
    **kwargs,
) -> scvi.model.SCVI:
    # copying the object so that we don't fuckup the source.  Could maybe add back an option to not do this to save
    # on memory
    adata_copy = adata.copy()

    # identify the HVGs and subset on those
    if isinstance(adata_copy, AnnData):
        scvi.data.poisson_gene_selection(adata_copy, layer=layers)
        mu.pp.filter_var(adata_copy, var="highly_variable")
    elif isinstance(adata_copy, MuData):
        if "rna" in adata_copy.mod:
            scvi.data.poisson_gene_selection(adata_copy[mudata_rna_modality], layer=layers[mudata_rna_modality])
            mu.pp.filter_var(adata_copy[mudata_rna_modality], var="highly_variable")

            # TODO: add columns for the isotype controls so that we can regress them out
            for _ in adata_copy[mudata_protein_modality].var_names[
                adata_copy[mudata_protein_modality].var_names.str.contains(isotype_pattern)
            ]:
                adata_copy.obs[_] = adata_copy[mudata_protein_modality].X[
                    :, adata_copy[mudata_protein_modality].var_names.get_loc(_)
                ]

    if layers and isinstance(adata_copy, AnnData):
        adata_copy.layers["X"] = adata_copy.X.copy()
        adata_copy.X = adata_copy.layers[layers].copy()
        adata_copy.layers[save_layer] = csr_matrix(adata_copy.X.copy())  # converts to CSR format, preserve counts
    if isinstance(adata_copy, MuData):
        for mod in adata_copy.mod:
            adata_copy[mod].layers[save_layer] = csr_matrix(
                adata_copy[mod].X.copy()
            )  # converts to CSR format, preserve counts

    # If we are working with RNA data, use scVI; for protein, use TotalVI
    match modality.lower():
        case "rna":
            scvi.model.SCVI.setup_anndata(
                adata_copy,
                layer=layers,
                batch_key=batch_key,
                size_factor_key=size_factor_key,
                categorical_covariate_keys=categorical_covariate_keys,
                continuous_covariate_keys=continuous_covariate_keys,
            )
            model = scvi.model.SCVI(adata_copy, dispersion=gene_dispersion, gene_likelihood=gene_likelihood)
        case "prot":
            # Two ways to set things up:
            # if an AnnData object, we need data in adata.obsm[protein_obsm_key]
            # and we need the feature names in adata.uns[protein_names_uns_key]
            if isinstance(adata, AnnData):
                if issparse(adata_copy.X):
                    adata_copy.obsm[protein_obsm_key] = adata_copy.X.toarray().copy()
                else:
                    adata_copy.obsm[protein_obsm_key] = adata_copy.X.copy()
                adata_copy.uns[protein_names_uns_key] = [f"{_}_x" for _ in adata_copy.var_names.to_list()]
                scvi.model.TOTALVI.setup_anndata(
                    adata=adata_copy,
                    layer=save_layer,
                    batch_key=batch_key,
                    protein_expression_obsm_key=protein_obsm_key,
                    # there is currently (2024-04-16) a bug in TotalVI that causes a
                    # "cannot reindex on an axis with duplicate labels" error.
                    # See: https://github.com/scverse/scvi-tools/issues/2627
                    # I, however, think this is because gene expression should stay in adata.X
                    # and prot data should be stuffed in adata.obsm[prot_key]
                    protein_names_uns_key=protein_names_uns_key,
                    size_factor_key=size_factor_key,
                    categorical_covariate_keys=categorical_covariate_keys,
                    continuous_covariate_keys=continuous_covariate_keys,
                )
            # things are much simpler if we have a MuData object
            elif isinstance(adata, MuData):
                if layers:
                    protein_layer = layers[mudata_protein_modality] if mudata_protein_modality in layers else None
                    rna_layer = layers[mudata_rna_modality] if mudata_rna_modality in layers else None
                else:
                    protein_layer = None
                    rna_layer = None
                isotype_keys = (
                    adata_copy[mudata_protein_modality].var_names[
                        adata_copy[mudata_protein_modality].var_names.str.contains(isotype_pattern)
                    ]
                    if use_isotype_controls
                    else None
                )
                if continuous_covariate_keys is None:
                    continuous_covariate_keys = []
                scvi.model.TOTALVI.setup_mudata(
                    mdata=adata_copy,
                    batch_key=batch_key,
                    rna_layer=rna_layer,
                    protein_layer=protein_layer,
                    continuous_covariate_keys=isotype_keys + continuous_covariate_keys,
                    modalities={"rna_layer": mudata_rna_modality, "protein_layer": mudata_protein_modality},
                    size_factor_key=size_factor_key,
                    categorical_covariate_keys=categorical_covariate_keys,
                )
            model = scvi.model.TOTALVI(
                adata_copy,
                protein_dispersion=protein_dispersion,
                gene_dispersion=gene_dispersion,
                gene_likelihood=gene_likelihood,
            )

    model.train(
        check_val_every_n_epoch=check_val_every_n_epoch,
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor=early_stopping_monitor,
        devices=num_devices,
        # **kwargs,
    )

    return model


def process_anndata(adata: AnnData, layers: str | list[str]) -> AnnData:
    if layers:
        layers = [layers] if isinstance(layers, str) else layers
        for _ in layers:
            if not is_integer_array(adata.layers[_]):
                msg = (
                    f"Both the scVI and TOTALVI models require raw, non-normalized counts and it looks like the {layers[_]} "
                    f"layer in the passed object are normalized"
                )
                raise ValueError(msg)
    elif not is_integer_array(adata.X):
        msg = (
            "the scVI and TOTALVI models require raw, non-normalized counts and it looks like the passed object has "
            "been normalized"
        )
        raise ValueError(msg)
    return adata


def process_mudata(
    mudata: MuData,
    layers: dict[str, str] | None = None,
) -> MuData:
    for mod in mudata.mod:
        if layers and (mod in layers and not is_integer_array(mudata[mod].layers[layers[mod]])):
            msg = (
                f"TOTALVI require raw, non-normalized counts and it looks like"
                f"layers[{mod}] in the {mod} modality is not all integers."
            )
            raise ValueError(msg)
        # not using a layer, so check X and if it isn't all integers, try pulling from raw if possible
        elif not is_integer_array(mudata[mod].X):
            if mudata[mod].raw and is_integer_array(mudata[mod].raw.X):
                mudata[mod].X = mudata[mod].raw[mudata[mod].obs_names, mudata[mod].var_names].X.copy()
                mudata[mod].X = mudata[mod].X.toarray() if issparse(mudata[mod].X) else mudata[mod].X
                # no? fuck it. error
            else:
                msg = (
                    f"TOTALVI requires raw, non-normalized counts and it looks a raw "
                    f"integer array for the {mod} modality could not be found in X nor the raw slot."
                )
                raise ValueError(msg)

    return mudata


def diff_expr_test(
    adata: AnnData | MuData,
    model: scvi.model.SCVI,
    groupby: str,
    current_group: str | None = None,
    reference_group: str | None = None,
    batch_correct: bool = False,
    remove_outliers: bool = False,
    **kwargs: Any,
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    match type(model):
        case scvi.model.SCVI:
            deg_df = model.differential_expression(
                groupby=groupby,
                group1=current_group,
                group2=reference_group,
                batch_correction=batch_correct,
                filter_outlier_cells=remove_outliers,
                **kwargs,
            )
        case scvi.model.TOTALVI:
            # need to change this when/if we add support for mudata objects.  It appears that TotalVI expects that
            # adata.X has gene expression while adata.obsm[protein_key] has the protein expression
            # if we do it correctly instead of how we're doing it here, we can get protein/gene correlations
            deg_df = model.differential_expression(
                groupby=groupby, group1=current_group, group2=reference_group, batch_correction=batch_correct, **kwargs
            )
        case _:
            msg = "Model is of an unrecognized type.  Currently only scVI and TotalVI models work."
            raise ValueError(msg)

    # does this even make sense since we rewrote to use the scvi-native 1 vs many?
    # if "group1" not in deg_df.columns or deg_df["group1"].nunique() <= 1:
    #     return deg_df
    if isinstance(adata, AnnData):
        return deg_df.loc[deg_df.index.intersection(adata.var_names), :]
    elif isinstance(adata, MuData):
        degs_dict = {i: deg_df.loc[deg_df.index.intersection(adata[i].var_names), :] for i in adata.mod}
    for _ in degs_dict:
        degs_dict[_] = degs_dict[_].loc[degs_dict[_].index.intersection(adata[_].var_names), :]
    return degs_dict


def is_integer_array(arr: npt.ArrayLike | csr_matrix) -> bool:
    """
    Test if an array is really all integers
    """
    if issparse(arr):
        # have to convert to `sparse.COO` here as `scipy.csr_matrix` either no longer or never did support numpy
        # calculations on it
        return not (np.mod(sparse.asCOO(arr), 1) != 0).any()
    else:
        return not (np.mod(arr, 1) != 0).any()
