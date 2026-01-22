############################################################################################################# 
# 使用scvi方法对批次进行矫正，画出批次矫正前后的umap图
############################################################################################################# 
def scvi(adata, outdir):
    import scanpy as sc
    import matplotlib.pyplot as plt
    import scvi
    
    # 查看数据是不是有批次效应
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata=adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata=adata, n_top_genes=3000, batch_key="Sample", subset=True)
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata=adata, n_pcs=50, log=True)
    sc.pp.neighbors(adata=adata, , n_pcs=30)
    sc.tl.umap(adata=adata)
    sc.pl.umap(adata, color=["Sample", "Class"], title=["Sample", "Class"], legend_loc=None) # 绘制批次校正前的UMAP图
    
    # scvi分析
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata=adata, layer="counts", batch_key="Sample")
    model = scvi.model.SCVI(adata=adata, n_layers=2, n_latent=30)
    model.train(enable_progress_bar=False)
    model.save(f"{outdir}/scvi_model", overwrite=True)
    
    # 查看批次效应去除后的效果
    model = scvi.model.SCVI.load(f"{outdir}/scvi_model", adata=adata)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata=adata, use_rep="X_scVI")
    sc.tl.umap(adata=adata)
    sc.pl.umap(adata, color=["Sample", "Class"], title=["Sample", "Class"], legend_loc=None) # 绘制批次校正后的UMAP图

    return(adata)
adata = scvi(adata=adata, outdir=dir)
adata.write_h5ad(f"{dir}/scvi.h5ad")
