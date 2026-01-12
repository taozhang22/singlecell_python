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
def multi_leiden_umapplot(adata, res_list, markers):
    import scanpy as sc

    for res in res_list:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
        sc.pl.dotplot(adata, markers, groupby=f"leiden_{res:.2f}", standard_scale="var", title=f"Resolution {res:.2f}")
    sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")

    return adata
adata = multi_leiden_umapplot(adata=adata, res_list=res_list, markers=markers)

############################################################################################################# 
# 计算每个簇的前十的差异表达基因
############################################################################################################# 
def deg(adata, my_res, outdir="result/base/"):
    import scanpy as sc
    
    sc.tl.rank_genes_groups(adata, groupby=f"leiden_{my_res:.2f}", method="wilcoxon")
    deg_df = sc.get.rank_genes_groups_df(adata, group=None)
    top10 = deg_df.groupby("group", as_index=False).head(10)
    top10.to_csv(f"{outdir}/top10_deg.csv", index=False)

    return adata
#####################################################使用案例#################################################################
# 暂停，根据情况更改项目参数
my_res = 0.3
adata = scell.base.deg(adata, my_res, outdir=dir)

# 将细胞进行注释
def celltype_annotation(
    adata, my_res, legend_fontsize=None, outdir="result/base/",
    ct_map={
        "0": "Epithelial",
        "1": "Plasma_B",
        "2": "T_NK",
        "3": "McDC",
        "4": "B_cell",
        "5": "Plasma_B",
        "6": "T_NK",
        "7": "Fibroblast",
        "8": "MAST",
        "9": "Enteric_glial",
        "10": "Endothelial"},
    markers={
        "B_cell": ["CD79A", "MS4A1", "BANK1"],
        "Endothelial": ["PLVAP", "VWF"],
        "Enteric_glial": ["PLP1", "CLU", "S100B"],
        "Epithelial": ["KRT18", "KRT8", "EPCAM"],
        "Fibroblast": ["COL1A1", "DCN", "COL3A1"],
        "MAST": ["MS4A2", "KIT", "TPSAB1"],
        "McDC": ["LYZ", "CD68", "CD14"],
        "Neutrophil": ["G0S2", "FCGR3B", "CSF3R"],
        "pDC": ["LILRA4", "IL3RA", "CLEC4C"],
        "Plasma_B": ["MZB1", "DERL3", "IGHG2"],
        "T_NK": ["CD3D", "CD3E", "KLRB1"]},
):
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
    fig = sc.pl.embedding(adata, basis="umap_after", color="Celltype", legend_loc="on data", show=False, return_fig=True)
    fig.savefig(f"{outdir}/pdf_umap_after_harmony_with_celltype.pdf", bbox_inches="tight")
    plt.show(); plt.close(fig)

    return adata



# 获取到含有完整count的adata对象
def adata_all_counts(adata,
                    file_adata_include_all_counts="result/base/qc.h5ad",
                    outdir_adata_include_all_counts="result/base/adata_all_counts.h5ad"):
    import scanpy as sc

    obs_df = adata.obs["Celltype"].copy()
    adata = sc.read_h5ad(file_adata_include_all_counts)
    adata.obs["Celltype"] = obs_df
    adata.write_h5ad(outdir_adata_include_all_counts)

#####################################################使用案例#################################################################
# 根据情况更改项目参数
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
adata = scell.base.multi_leiden(adata, res_list=res_list, markers=markers)
adata.write_h5ad(f"{dir}/adata_leiden.h5ad")

# 暂停，根据情况更改项目参数
my_res = 0.3
adata = scell.base.deg(adata, my_res, outdir=dir)
