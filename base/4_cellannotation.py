# 参数，根据实际情况进行修改---------------------------------------
res_list = np.arange(0.1, 1.1, 0.1)
markers = {
    "B_cell": ["CD79A", "MS4A1", "BANK1"],
    "Endothelial": ["PLVAP", "VWF", "CLDN5"],
    "Enteric_glial": ["PLP1", "CLU", "S100B"],
    "Epithelial": ["KRT18", "KRT8", "EPCAM"],
    "Fibroblast": ["COL1A1", "DCN", "COL3A1"],
    "MAST": ["MS4A2", "KIT", "TPSAB1"],
    "McDC": ["LYZ", "CD68", "CD14"],
    "Neutrophil": ["G0S2", "FCGR3B", "CSF3R"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C"],
    "Pericyte": ["MCAM", "PDGFRB", "CSPG4"],
    "Plasma_B": ["MZB1", "DERL3", "IGHG2"],
    "SMC": ["CNN1", "TAGLN", "DES"],  # Smooth_muscle_cell
    "T_NK": ["CD3D", "CD3E", "KLRB1"],
}

for res in res_list:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
    sc.pl.dotplot(adata, markers, groupby=f"leiden_{res:.2f}", standard_scale="var", title=f"Resolution {res:.2f}")
sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")


# 参数，根据实际情况进行修改---------------------------------------------------
my_res = 0.3

sc.tl.rank_genes_groups(adata, groupby=f"leiden_{my_res:.2f}", method="wilcoxon", pts=True)
deg = sc.get.rank_genes_groups_df(adata, group=None)
deg = deg.query("pvals_adj < 0.05 and logfoldchanges > 0.585 and pct_nz_group >= 0.25")
top20 = deg.groupby("group", as_index=False).head(20)
top20.to_csv(f"{dir}/top_deg.csv", index=False)


# 参数，根据实际情况进行修改
ct_map = {
    "0": "Epithelial",
    "1": "T_NK",
    "2": "T_NK",
    "3": "McDC",
    "4": "Epithelial",
    "5": "B_cell",
    "6": "Fibroblast",
    "7": "Plasma_B",
    "8": "Fibroblast",
    "9": "Endothelial",
    "10": "Plasma_B",
    "11": "Plasma_B",
    "12": "SMC",
    "13": "Enteric_glial",
    "14": "MAST",
    "15": "pDC",
    "16": "Endothelial",
    "17": "Neutrophil",
}

adata.obs["Celltype"] = adata.obs[f"leiden_{my_res:.2f}"].map(ct_map).astype("category") # 增加Celltype列
cats_sorted = sorted(adata.obs["Celltype"].cat.categories, key=str.lower) # 获取Celltype的顺序
adata.obs["Celltype"] = adata.obs["Celltype"].cat.reorder_categories(cats_sorted, ordered=True) # 将Celltype类别的顺序设置为不区分大小写，按照字母顺序排序
celltypes = set(adata.obs["Celltype"].unique()) # 获取被注释到的细胞
filtered_genes = {ct: mk for ct, mk in markers.items() if ct in celltypes} # 保留被注释到的细胞和基因的字典

# 绘制注释后的气泡图
sc.pl.dotplot(adata, var_names=filtered_genes, groupby="Celltype", standard_scale="var", show=False)
plt.savefig(f"{dir}/dotplot_after_celltype_annotation.pdf", bbox_inches="tight")
plt.show(); plt.close()

# 绘制注释后的umap图
fig = sc.pl.umap(adata, color="Celltype", legend_loc="on data", show=False, return_fig=True, legend_fontsize=legend_fontsize)
fig.savefig(f"{dir}/umap_after_celltype_annotation.pdf", bbox_inches="tight")
plt.show(); plt.close(fig)

adata.write_h5ad(f"{dir}/celltype_annotated.h5ad")


# 制作完整count的adata----------------------------------------------
obs_df = adata.obs["Celltype"].copy()
adata = sc.read_h5ad(f"{dir}/qc.h5ad")
adata.obs["Celltype"] = obs_df.reindex(adata.obs_names)
adata.write_h5ad(f"{dir}/adata.h5ad")
