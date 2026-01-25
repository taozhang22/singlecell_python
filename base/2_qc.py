############################################################################################################# 
# 对数据进行质量控制
############################################################################################################# 
min_genes=100; min_cells=3; pct_counts_mt=20

# 计算质量控制指标
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)

# 绘制质控图片
fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
plt.tight_layout(); plt.show(); plt.close(fig)
    
fig, axes = plt.subplots(5, 1, figsize=(0.3 * adata.obs["Sample"].nunique(), 15))
for i, key in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]):
    sc.pl.violin(adata, keys=key, groupby="Sample", jitter=0.4, rotation=45, show=False, ax=axes[i])
plt.tight_layout(); plt.show(); plt.close(fig)

# 质控
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
adata = adata[adata.obs["pct_counts_mt"] < pct_counts_mt].copy()

# 去除双细胞
sc.pp.scrublet(adata, batch_key="Sample")
adata = adata[~adata.obs["predicted_doublet"]].copy()
adata.write_h5ad(f"{dir}/qc.h5ad")
