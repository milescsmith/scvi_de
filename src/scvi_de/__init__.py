import warnings
from typing import Literal

import anndata as ad
import pandas as pd
import scvi
from scipy.sparse import csr_matrix, issparse


# TODO: pass more options on
def scvi_de(
    adata: ad.AnnData | None = None,
    model: scvi.model.SCVI | scvi.model.TOTALVI | None = None,
    groupby: str = "leiden",
    layer: str | None = None,
    compare_group: str | None = None,
    reference_group: str = "rest",
    batch_correct: bool = False,
    batch_key: str = "batch",
    dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
    gene_likelihood: Literal["nb", "zinb", "poisson"] = "nb",
    remove_outliers: bool = False,
    inplace: bool = False,
    return_df: bool = False,
    return_model: bool = False,
    lfc_use: Literal["lfc_mean", "lfc_median", "lfc_max", "lfc_min"] = "lfc_mean",
    **kwargs,
) -> dict[str, pd.DataFrame | ad.AnnData | scvi.model.SCVI | scvi.model.TOTALVI] | None:
    """Perform differential expression analysis using the scVI suite of tools

    Parameters
    -------
    adata : anndata.AnnData
        Object to model and perform differential expression on. Either an AnnData or model object need to be passed;
        if both are passed, AnnData will be ignored.
    model : :class:`~scvi.model.SCVI`, :class:`~scvi.model.TOTALVI`
        A previously trained model.  If both an AnnData and model object are passed, the AnnData object will be ignored
        and the model will be used.
    groupby : str, default="leiden"
        Column in `obs` to use when performing comparisons.
    layer : str
        Optional. Key for `adata.layers` where raw count data is stored. If not provided, the counts in `adata.X` will
        be used
    compare_group : str
        Optional. Category in `adata.obs[groupby]` to use as the query for comparison of two groups. Positive fold
        changes indicate higher expression in this group.
    reference_group : str
        Optional. Category in `adata.obs[groupby]` to use as a reference for comparison of two groups. Negative
        fold changes indicate higher expression in this group.
    batch_correct : bool, default=False
        Optional. Should batch effects be corrected for in DE inference?
    batch_key : str, default="batch"
        Optional: Column in `adata.obs` to use to indicate batches.
    dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"], default="gene"
        One of the following:

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    gene_likelihood : Literal["nb", "zinb", "poisson"], default="nb"
        One of the following:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution

    remove_outliers : bool, default=False
        Whether to filter outlier cells with `:meth:~scvi.model.base.DifferentialComputation.filter_outlier_cells()`
    inplace : bool, default=False
        Should any changes be performed directly on the provided AnnData object. If not, a copy is made.
    return_df : bool, default=False
        Should a dataframe of results be returned?
    return_model : bool, default=False
        Should the model be returned?
    lfc_use : Literal["lfc_mean", "lfc_median", "lfc_max", "lfc_min"], default="lfc_mean"
        Which log fold change to use.

    Return
    ------
        dict of {str: pd.DataFrame | ad.AnnData | scvi.model.SCVI | scvi.model.TOTALVI}
    """
    if adata and model:
        msg = (
            "Passing and Anndata object and scVI model is confusing. "
            "If you wish to regenerate the model, pass just an Anndata object. "
            "For now, using the previous model in the interests of time..."
        )
        warnings.warn(msg, UserWarning, stacklevel=2)
    elif adata:
        model = create_model(
            adata=adata,
            layer=layer,
            batch_key=batch_key,
            inplace=inplace,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            **kwargs,
        )
    else:
        msg = "Either an AnnData or scVI model needed and both are missing."
        raise ValueError(msg)

    adata_copy = model.adata
    latent = model.get_latent_representation()
    adata_copy.obsm["X_scVI"] = latent
    if compare_group:
        de_change = one_vs_others(
            adata=adata_copy,
            model=model,
            groupby=groupby,
            current_group=compare_group,
            reference_group=reference_group,
            batch_correct=batch_correct,
            remove_outliers=remove_outliers,
        )
    else:
        de_change = {
            i: one_vs_others(
                adata=adata_copy,
                model=model,
                groupby=groupby,
                current_group=i,
                batch_correct=batch_correct,
                remove_outliers=remove_outliers,
            )
            for i in adata_copy.obs[groupby].unique()
        }

    names = {
        i: (de_change[i]["proba_de"] * de_change[i][lfc_use] * de_change[i]["non_zeros_proportion1"])
        .sort_values(ascending=False)
        .index.to_list()
        for i in sorted(de_change)
    }

    scores_rec = pd.DataFrame(
        {
            i: (de_change[i]["proba_de"] * de_change[i][lfc_use] * de_change[i]["non_zeros_proportion1"])
            .sort_values(ascending=False)
            .to_numpy()
            for i in sorted(de_change)
        }
    ).to_records(index=False)

    pvals_rec = pd.DataFrame(
        {i: de_change[i].loc[names[i], "proba_not_de"].to_numpy() for i in sorted(de_change)}
    ).to_records(index=False)

    lfc_rec = pd.DataFrame({i: de_change[i].loc[names[i], lfc_use].to_numpy() for i in sorted(de_change)}).to_records(
        index=False
    )

    de_dict = {
        "params": {
            "groupby": "leiden_wnn_0.9",
            "reference": "rest",
            "method": "scvi_de",
            "use_raw": True,
            "layer": None,
            "lfc_used": lfc_use,
            "corr_method": "benjamini-hochberg",
        },
        "names": pd.DataFrame(names).to_records(index=False),
        "scores": scores_rec,
        "pvals": pvals_rec,
        "pvals_adj": pvals_rec,
        "logfoldchanges": lfc_rec,
    }

    return_dict = {}
    if inplace:
        adata.uns["rank_genes_groups"] = de_dict
    else:
        adata_copy.uns["rank_genes_groups"] = de_dict
        return_dict["adata"] = adata_copy

    if return_df:
        return_dict["deg_df"] = de_change
    if return_model:
        return_dict["model"] = model
    return return_dict


def create_model(
    adata: ad.AnnData,
    layer: str | None = None,
    model_type: Literal["scvi", "totalvi"] = "scvi",  # add more options and we validate input arguments
    batch_key: str = "str",
    inplace: bool = False,
    dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
    gene_likelihood: Literal["nb", "zinb", "poisson"] = "nb",
    save_layer: str = "counts",
    check_val_every_n_epoch: int = 1,
    max_epochs: int = 400,
    early_stopping: bool = True,
    early_stopping_patience: int = 20,
    early_stopping_monitor: Literal[
        "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
    ] = "elbo_validation",
    **kwargs,
) -> scvi.model.SCVI:
    adata_copy = adata if inplace else adata.copy()
    if layer:
        if ((adata_copy.layers[layer] % 0) != 0).any():
            msg = f"The {model_type} model requires raw, non-normalized counts and it looks like the {layer} layer in the passed object are normalized"
            raise ValueError(msg)
    elif issparse(adata_copy.X):
        if ((adata_copy.X.toarray() % 0) != 0).any():
            msg = f"The {model_type} model requires raw, non-normalized counts and it looks like the passed object has been normalized"
            raise ValueError(msg)
    elif ((adata_copy.X % 0) != 0).any():
        msg = f"The {model_type} model requires raw, non-normalized counts and it looks like the passed object has been normalized"
        raise ValueError(msg)

    scvi.data.poisson_gene_selection(adata_copy, layer=layer)

    adata_copy = adata_copy[:, adata_copy.var["highly_variable"]]  # focus on selected genes

    adata_copy.layers[save_layer] = csr_matrix(adata_copy.X.copy())  # converts to CSR format, preserve counts

    match model_type.lower():
        case "scvi":
            scvi.model.SCVI.setup_anndata(adata_copy, layer=layer, batch_key=batch_key)
            model = scvi.model.SCVI(adata_copy, dispersion=dispersion, gene_likelihood=gene_likelihood)
        case "totalvi":
            # assuming that we're working with the protein portion of a MuData object
            adata_copy.obsm["prot"] = adata_copy.X
            adata.uns["prot_names"] = adata.var_names.to_list()
            scvi.model.TOTALVI.setup_anndata(
                adata=adata_copy,
                layer=save_layer,
                batch_key=batch_key,
                protein_expression_obsm_key="prot",
                protein_names_uns_key="prot_names",
            )
            model = scvi.model.TOTALVI(adata_copy, dispersion=dispersion, gene_likelihood=gene_likelihood)

    model.train(
        check_val_every_n_epoch=check_val_every_n_epoch,
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor=early_stopping_monitor,
        **kwargs,
    )

    return model


def one_vs_others(
    adata: ad.AnnData,
    model: scvi.model.SCVI,
    groupby: str,
    current_group: str,
    reference_group: str = "rest",
    batch_correct: bool = False,
    remove_outliers: bool = False,
):
    cell_idx1 = adata.obs[groupby] == current_group

    if reference_group == "rest":
        cell_idx2 = adata.obs[groupby] != current_group
    else:
        cell_idx2 = adata.obs[groupby] == reference_group

    return model.differential_expression(
        idx1=cell_idx1, idx2=cell_idx2, batch_correction=batch_correct, filter_outlier_cells=remove_outliers
    )
