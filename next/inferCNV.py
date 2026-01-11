############################################################################################################# 
# 准备数据
############################################################################################################# 
dir = "result/inferCNV"
os.makedirs(dir, exist_ok=True)

tumor = ['Tumor', 'Border'] # 需要根据情况修改参数
adata = sc.read_h5ad("result/base/adata.h5ad")
adata = adata[adata.obs["Celltype"] == "Epithelial"].copy()
adata_scvi = adata[adata.obs['Class'].isin(tumor)].copy()

############################################################################################################# 
# 将数据进行分簇，然后绘制出umap图
############################################################################################################# 
def multi_leiden_umapplot(adata, res_list, dotplot):
    import scanpy as sc

    for res in res_list:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
    sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")

    return adata
