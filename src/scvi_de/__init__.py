from typing import Literal

import anndata as ad
import pandas as pd
import scvi
from scipy.sparse import csr_matrix


# TODO: pass more options on
def scvi_de(
    adata: ad.AnnData,
    groupby: str,
    compare_group: str | None = None,
    reference_group: str = "rest",
    batch_correct: bool = False,
    batch_key: str = "batch",
    remove_outliers: bool = False,
    # plot: bool = False,
    inplace: bool = False,
    return_df: bool = False,
    return_model: bool = False,
    lfc_use: Literal["lfc_mean", "lfc_median", "lfc_max", "lfc_min"] = "lfc_mean",
) -> dict[str, pd.DataFrame | ad.AnnData | scvi.model.SCVI] | None:
    adata_copy = adata.copy()

    scvi.data.poisson_gene_selection(adata_copy)

    adata_copy = adata_copy[:, adata_copy.var["highly_variable"]]  # focus on selected genes

    adata_copy.layers["counts"] = csr_matrix(adata_copy.X.copy())  # converts to CSR format, preserve counts

    scvi.model.SCVI.setup_anndata(adata_copy, layer="counts", batch_key=batch_key)

    model = scvi.model.SCVI(adata_copy, gene_likelihood="nb")
    model.train(
        check_val_every_n_epoch=1,
        max_epochs=400,
        early_stopping=True,
        early_stopping_patience=20,
        early_stopping_monitor="elbo_validation",
    )

    # I'll renable this when I feel like adding matplotlib to the list of dependencies
    # if plot:
    #     # Ensure convergence
    #     train_test_results = model.history["elbo_train"]
    #     train_test_results["elbo_validation"] = model.history["elbo_validation"]
    #     train_test_results.iloc[10:].plot(logy=True)  # exclude first 10 epochs
    #     plt.show()

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

    if inplace:
        adata.uns["rank_genes_groups"] = de_dict
        return None
    else:
        adata_copy.uns["rank_genes_groups"] = de_dict
        return_dict = {"adata": adata_copy}

        if return_df:
            return_dict["deg_df"] = de_change
        if return_model:
            return_dict["model"] = model
        return return_dict


def one_vs_others(
    adata: ad.AnnData,
    model,
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
