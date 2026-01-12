# https://infercnvpy.readthedocs.io/en/latest/infercnv.html
############################################################################################################# 
# 设置参数
#############################################################################################################
dir = "result/inferCNV"
tumor = ['Tumor', 'Border']
normal = ["Normal"]
adata_input = "result/base/adata.h5ad"
my_res = 0.4

############################################################################################################# 
# 准备
#############################################################################################################
def prepare_infercnv(
    dir="result/inferCNV",
    adata_input="result/base/adata.h5ad",
    tumor=['Tumor', 'Border']):
    import os
    import scanpy as sc
    os.makedirs(dir, exist_ok=True)

    adata = sc.read_h5ad(adata_input)
    adata = adata[adata.obs["Celltype"] == "Epithelial"].copy()
    adata_scvi = adata[adata.obs['Class'].isin(tumor)].copy()
adata_scvi = prepare_infercnv(dir=dir, adata_input=adata_input, tumor=tumor)

############################################################################################################# 
# 去除批次效应
#############################################################################################################
# 此处引用base文件夹下3_batch_correction文件的内容

############################################################################################################# 
# 将数据进行分簇，然后绘制出umap图
############################################################################################################# 
def multi_leiden_umapplot(adata, res_list):
    import scanpy as sc

    for res in res_list:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
    sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")

    return adata
adata = multi_leiden_umapplot(adata=adata_scvi, res_list=np.arange(0.1, 1.1, 0.1))
