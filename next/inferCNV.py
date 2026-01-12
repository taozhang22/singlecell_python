# https://infercnvpy.readthedocs.io/en/latest/infercnv.html
############################################################################################################# 
# 设置参数
#############################################################################################################
# 全局参数，根据情况修改
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

############################################################################################################# 
# 制作cnv分析需要的label标签，读入染色体位置文件，后进行分析，画出染色体cnv图
############################################################################################################# 
adata_tumor = adata[adata.obs["Class"].isin(tumor)].copy()
adata_tumor.obs["label"] = adata_scvi.obs[f"leiden_{my_res:.2f}"].reindex(adata_tumor.obs_names)
adata_normal = adata[adata.obs["Class"].isin(normal)].copy()
adata_normal.obs["label"] = "Control"
adata = sc.concat([adata_normal, adata_tumor])

var_df = pd.read_csv("../../resource/inferCNV/hg38_gencode_v27.txt", sep="\t", index_col=0, header=None, names=["chromosome", "start", "end"])
adata.var = var_df.join(adata.var, how="right")

cnv.tl.infercnv(adata, reference_key="label", reference_cat=["Control"], window_size=250)
ax = cnv.pl.chromosome_heatmap(adata, groupby="label", dendrogram=False)
ax.figure.savefig(f"{dir}/chromosome_heatmap_label.png")
plt.show(); plt.close()

############################################################################################################# 
# 计算每个组的平均cnv分数
############################################################################################################# 
cnv.tl.cnv_score(adata, groupby="label")
display(adata.obs[["label", "cnv_score"]].drop_duplicates())

############################################################################################################# 
# 计算每个细胞的平均cnv分数，并且绘制箱线图
############################################################################################################# 
x = adata.obsm["X_cnv"].toarray()
x = StandardScaler().fit_transform(x)
x = MinMaxScaler(feature_range=(-1, 1)).fit_transform(x)
adata.obs["cnv_score_cell"] = (x * x).sum(axis=1)

order = (adata.obs.groupby("label")["cnv_score_cell"].median().sort_values(ascending=False).index.tolist())
order.remove("Control")
order.append("Control")
sns.boxplot(data=adata.obs, x="label", y="cnv_score_cell", order=order, showfliers=False)
