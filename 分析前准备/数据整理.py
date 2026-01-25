# Environment: scenv
import os
import pandas as pd
import datatable as dt
import matplotlib.pyplot as plt
import scanpy as sc
os.chdir("/home/yang/research/bioinformation/singlecell/database") # 根据情况修改工作路径

# 参数，根据实际情况修改
dir = "GSE132465"
file_count = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
file_meta1 = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"

# 读取表达矩阵，形成adata
count = dt.fread(f"{dir}/{file_count}", sep="\t")
adata = count[:, 1:].to_pandas() # matrix
adata.index = count[:, 0].to_list()[0] # index
adata = adata.T

adata=sc.AnnData(adata)
adata.var_names_make_unique()

# 读取注释信息
meta1 = pd.read_csv(f"{dir}/{file_meta1}", index_col=0, sep="\t")
meta1 = meta1[["Sample", "Class", "Cell_type", "Cell_subtype"]]

meta2 = pd.read_csv(f"{dir}/{dir}.txt", sep="\t", index_col=0)
meta2 = meta2[["MSI"]]

meta = meta1.join(meta2, on="Sample", how="inner")
display(meta)

# 形成adata文件
adata.obs = meta.copy()
adata.write_h5ad(f"{dir}/python/{dir}.h5ad")
adata.obs["Sample"].value_counts()

# 质量控制
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)

for Sample in adata.obs["Sample"].astype(str).unique().tolist():
    fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))
    sc.pl.scatter(adata[adata.obs['Sample'] == Sample], x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
    sc.pl.scatter(adata[adata.obs['Sample'] == Sample], x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
    plt.tight_layout(); plt.show(); plt.close(fig)
    
fig, axes = plt.subplots(5, 1, figsize=(0.3 * adata.obs["Sample"].nunique(), 15))
for i, key in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]):
    sc.pl.violin(adata, keys=key, groupby="Sample", jitter=0.4, rotation=45, show=False, ax=axes[i])
plt.tight_layout(); plt.show(); plt.close(fig)

adata.obs["Sample"].value_counts()

