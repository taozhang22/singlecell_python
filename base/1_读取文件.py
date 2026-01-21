# 使用说明本：本代码的细胞注释之前的部分要修改的地方很少，细胞注释阶段有多个参数需要修改，请仔细阅读注释部分

############################################################################################################# 
# 读取本研究需要的adata文件
#############################################################################################################
# 方法一：
def read_files(filename, pattern):
    from pathlib import Path
    import scanpy as sc

    files = list(Path(filename).rglob(pattern))
    for file in files: print(file)
    adata = sc.concat([sc.read_h5ad(f) for f in files])

    return adata
adata = read_files(filename="/home/students/zhangtao/research/singlecell/database", pattern="*(韩国).h5ad") # 需要根据实际的工作路径和文件名称修改

# 方法二
patterns = ["GSE132465", "GSE144735"] # 根据实际想要纳入的数据集名称修改
filenames = [f"../../database/{pattern}/python/{pattern}.h5ad" for pattern in patterns]
adata = sc.concat([sc.read_h5ad(f) for f in filenames])
