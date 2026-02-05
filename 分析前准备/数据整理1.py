%load_ext autotime
import warnings
warnings.filterwarnings("ignore")

import os
from pathlib import Path
import pandas as pd
import datatable as dt
import matplotlib.pyplot as plt
import scanpy as sc
os.chdir("/home/yang/research/bioinformation/singlecell/database") # 根据实际工作路径进行修改

# 以下参数根据实际情况进行修改
dir = "GSE132465"
file_count = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
file_meta1 = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"

# 制作adata
count = dt.fread(f"{dir}/{file_count}", sep="\t")
adata = count[:, 1:].to_pandas()
adata.index = count[:, 0].to_list()[0]
adata = adata.T

adata=sc.AnnData(adata)
adata.var_names_make_unique()

# 制作注释信息文件
meta1 = pd.read_csv(f"{dir}/{file_meta1}", index_col=0, sep="\t")
meta1 = meta1[["Sample", "Class"]]

meta2 = pd.read_csv(f"{dir}/{dir}.txt", sep="\t", index_col=0)
meta2 = meta2[["MSI"]]

meta = meta1.join(meta2, on="Sample", how="inner")
adata.obs = meta.copy()

###################################################################################
p = Path("GSE231559/data")
seen = set()
for f in p.iterdir():
    name = f.name
    if "_barcodes" in name:
        seen.add(name.split("barcodes", 1)[0])
    elif "_features" in name:
        seen.add(name.split("features", 1)[0])
    elif "_matrix" in name:
        seen.add(name.split("matrix", 1)[0])
prefixes = sorted(seen)

adatas = []
for prefix in prefixes:
    adata = sc.read_10x_mtx(path="GSE231559/data", var_names="gene_symbols", cache=True, prefix=prefix)
    adata.obs["Sample"] = prefix.split("_")[0]
    adatas.append(adata)
adata = sc.concat(adatas)

###################################################################################################
dir = "GSE178318"
mtx = sc.read_mtx(f"{dir}/GSE178318_matrix.mtx.gz").T
barcodes = pd.read_csv(f"{dir}/GSE178318_barcodes.tsv.gz", header=None, sep="\t", index_col=0)
barcodes.index.name = None
genes = pd.read_csv(f"{dir}/GSE178318_genes.tsv.gz", header=None, sep="\t", index_col=1)
genes.index.name = None
genes.columns = ["gene_id"]

adata = sc.AnnData(X=mtx.X, obs=barcodes, var=genes)
adata.obs["Sample"] = ["_".join(x.split("_")[1:3]) for x in adata.obs_names]
adata.obs.index.name = None
adata.write_h5ad(f"{dir}/python/{dir}.h5ad")























# 质量控制图
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)

for Sample in adata.obs["Sample"].astype(str).unique().tolist():
    fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))
    adat = adata[adata.obs['Sample'] == Sample].copy()
    sc.pl.scatter(adat, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
    sc.pl.scatter(adat, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
    plt.tight_layout(); plt.show(); plt.close(fig)

fig, axes = plt.subplots(5, 1, figsize=(0.3 * adata.obs["Sample"].nunique(), 15))
for i, key in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]):
    sc.pl.violin(adata, keys=key, groupby="Sample", jitter=0.4, rotation=45, show=False, ax=axes[i])
plt.tight_layout(); plt.show(); plt.close(fig)
print(adata.obs["Sample"].value_counts())

# 删除细胞数少于1000的样本
ncell = adata.obs["Sample"].astype(str).value_counts()
adata.obs["ncell"] = adata.obs["Sample"].astype(str).map(ncell)
adata = adata[adata.obs["ncell"] > 1000].copy()
adata.write_h5ad(f"{dir}/python/{dir}.h5ad")
