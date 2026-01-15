############################################################################################################# 
# 参数
############################################################################################################# 
# 细胞注释阶段用到的全局变量，根据研究的实际情况进行修改
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

############################################################################################################# 
# 算出umap图中分为多少个簇，并且绘制气泡图和umap图
############################################################################################################# 
def multi_leiden_dotplot_umapplot(adata, res_list, markers):
    import scanpy as sc

    for res in res_list:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
        sc.pl.dotplot(adata, markers, groupby=f"leiden_{res:.2f}", standard_scale="var", title=f"Resolution {res:.2f}")
    sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")

    return adata
adata = multi_leiden_dotplot_umapplot(adata=adata, res_list=res_list, markers=markers)

############################################################################################################# 
# 参数
############################################################################################################# 
# 细胞注释阶段用到的全局变量，根据研究的实际情况进行修改
my_res = 0.3

############################################################################################################# 
# 计算每个簇的前十的差异表达基因
############################################################################################################# 
def deg(adata, my_res, outdir):
    import scanpy as sc
    
    sc.tl.rank_genes_groups(adata, groupby=f"leiden_{my_res:.2f}", method="wilcoxon")
    deg_df = sc.get.rank_genes_groups_df(adata, group=None)
    top10 = deg_df.groupby("group", as_index=False).head(10)
    top10.to_csv(f"{outdir}/top10_deg.csv", index=False)

    return adata
adata = deg(adata=adata, my_res=my_res, outdir=dir)

############################################################################################################# 
# 参数
#############################################################################################################
# 细胞注释阶段用到的全局变量，根据研究的实际情况进行修改
# 强烈注意，每跑一次都要修改
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

############################################################################################################# 
# 注释，画出注释完毕的气泡图和umap图
############################################################################################################# 
# 将细胞进行注释
def celltype_annotation(adata, my_res, legend_fontsize, outdir, ct_map, markers):
    import scanpy as sc
    import matplotlib.pyplot as plt

    # 将Celltype（属于categories）的类别按照字母顺序排序，不区分大小写
    adata.obs["Celltype"] = adata.obs[f"leiden_{my_res:.2f}"].map(ct_map).astype("category") # 增加Celltype列
    cats_sorted = sorted(adata.obs["Celltype"].cat.categories, key=str.lower) # 获取Celltype的顺序
    adata.obs["Celltype"] = adata.obs["Celltype"].cat.reorder_categories(cats_sorted, ordered=True) # 将Celltype类别的顺序设置为不区分大小写，按照字母顺序排序
    celltypes = set(adata.obs["Celltype"].unique()) # 获取被注释到的细胞
    filtered_genes = {ct: mk for ct, mk in markers.items() if ct in celltypes} # 保留被注释到的细胞和基因的字典

    # 绘制注释后的气泡图
    sc.pl.dotplot(adata, var_names=filtered_genes, groupby="Celltype", standard_scale="var", show=False)
    plt.savefig(f"{outdir}/dotplot_after_annotation.pdf", bbox_inches="tight")
    plt.show(); plt.close()

    # 绘制注释后的umap图
    fig = sc.pl.umap(adata, color="Celltype", legend_loc="on data", show=False, return_fig=True, legend_fontsize=legend_fontsize)
    fig.savefig(f"{outdir}/pdf_umap_after_harmony_with_celltype.pdf", bbox_inches="tight")
    plt.show(); plt.close(fig)

    return adata
adata = celltype_annotation(adata=adata, my_res=my_res, legend_fontsize=None, outdir=dir, ct_map=ct_map, markers=markers)
adata.write_h5ad(f"{dir}/celltype_annotated.h5ad")

############################################################################################################# 
# 注释，画出注释完毕的气泡图和umap图
############################################################################################################# 
def adata_all_counts(adata, file_adata_include_all_counts, outdir_adata_include_all_counts):
    import scanpy as sc

    obs_df = adata.obs["Celltype"].copy()
    adata = sc.read_h5ad(file_adata_include_all_counts)
    adata.obs["Celltype"] = obs_df
    adata.write_h5ad(outdir_adata_include_all_counts)
adata_all_counts(adata=adata, file_adata_include_all_counts=f"{dir}/qc.h5ad", outdir_adata_include_all_counts=f"{dir}/adata.h5ad")
