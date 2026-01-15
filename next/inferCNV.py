# https://infercnvpy.readthedocs.io/en/latest/infercnv.html
############################################################################################################# 
# 设置参数
#############################################################################################################
# 全局参数，根据情况修改
import infercnvpy as cnv
dir = "result/inferCNV"
os.makedirs(dir, exist_ok=True)
tumor = ['Tumor', 'Border']
normal = ["Normal"]
adata_input = "result/base/adata.h5ad"
my_res = 0.4

############################################################################################################# 
# 准备，获取文件
#############################################################################################################
adata = sc.read_h5ad(adata_input)
adata = adata[adata.obs["Celltype"] == "Epithelial"].copy()
adata_scvi = adata[adata.obs['Class'].isin(tumor)].copy()

############################################################################################################# 
# 去除批次效应，引用的base文件夹的去除批次效应文件
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
    sc.pp.neighbors(adata=adata)
    sc.tl.umap(adata=adata)
    sc.pl.umap(adata, color=["Sample", "Class"], title=["Sample", "Class"], legend_loc=None) # 绘制批次校正前的UMAP图
    
    # scvi分析
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata=adata, layer="counts", batch_key="Sample")
    model = scvi.model.SCVI(adata=adata, n_layers=2, n_latent=30)
    model.train()
    model.save(f"{outdir}/scvi_model", overwrite=True)
    
    # 查看批次效应去除后的效果
    model = scvi.model.SCVI.load(f"{outdir}/scvi_model", adata=adata)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata=adata, use_rep="X_scVI")
    sc.tl.umap(adata=adata)
    sc.pl.umap(adata, color=["Sample", "Class"], title=["Sample", "Class"], legend_loc=None) # 绘制批次校正后的UMAP图

    return(adata)
adata_scvi = scvi(adata=adata_scvi, outdir=dir)

############################################################################################################# 
# 将数据进行分簇，然后绘制出umap图
############################################################################################################# 
def multi_leiden_umapplot(adata, res_list):
    import scanpy as sc

    for res in res_list:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
    sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")

    return adata
adata_scvi = multi_leiden_umapplot(adata=adata_scvi, res_list=np.arange(0.1, 1.1, 0.1))

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
ax_dict = cnv.pl.chromosome_heatmap(adata, groupby="label", dendrogram=True, show=False)
fig = next(iter(ax_dict.values())).figure
fig.savefig(f"{dir}/chromosome_heatmap_label.pdf", bbox_inches="tight")

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
adata.write_h5ad(f"{dir}/infercnv.h5ad")

obs_df = adata.obs[["label", "cnv_score_cell"]].copy()
obs_df["label"] = obs_df["label"].astype(str)
orders = (obs_df.groupby("label")["cnv_score_cell"].median().sort_values(ascending=False).index.tolist())
orders.remove("Control")
orders.append("Control")

sns.boxplot(data=obs_df, x="label", y="cnv_score_cell", order=orders, showfliers=False)
plt.xlabel("Group")
plt.ylabel("CNV Score")

############################################################################################################# 
# 计算各个组的差异分析的p值
############################################################################################################# 
from scipy.stats import mannwhitneyu
U, p = mannwhitneyu(
    x = obs_df[obs_df["label"] != "0"]["cnv_score_cell"],
    y = obs_df[obs_df["label"] == "Control"]["cnv_score_cell"],
    alternative="two-sided")
print(f"Mann-Whitney U test p-value: {p}")

U, p = mannwhitneyu(
    x = obs_df[obs_df["label"] != "1"]["cnv_score_cell"],
    y = obs_df[obs_df["label"] == "Control"]["cnv_score_cell"],
    alternative="two-sided")
print(f"Mann-Whitney U test p-value: {p}")

U, p = mannwhitneyu(
    x = obs_df[obs_df["label"] != "2"]["cnv_score_cell"],
    y = obs_df[obs_df["label"] == "Control"]["cnv_score_cell"],
    alternative="two-sided")
print(f"Mann-Whitney U test p-value: {p}")

U, p = mannwhitneyu(
    x = obs_df[obs_df["label"] != "3"]["cnv_score_cell"],
    y = obs_df[obs_df["label"] == "Control"]["cnv_score_cell"],
    alternative="two-sided")
print(f"Mann-Whitney U test p-value: {p}")

U, p = mannwhitneyu(
    x = obs_df[obs_df["label"] != "4"]["cnv_score_cell"],
    y = obs_df[obs_df["label"] == "Control"]["cnv_score_cell"],
    alternative="two-sided")
print(f"Mann-Whitney U test p-value: {p}")
