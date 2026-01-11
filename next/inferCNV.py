
def multi_leiden_dotplot_umapplot(adata, res_list, dotplot):
    import scanpy as sc

    for res in res_list:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
    sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")

    return adata
  
