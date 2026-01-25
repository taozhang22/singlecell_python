# Environment: scenv
import os
import pandas as pd
import scanpy as sc
import datatable as dt
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


