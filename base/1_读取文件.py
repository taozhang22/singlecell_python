############################################################################################################# 
# 设置整个代码多次需要的参数，放在设置环境的单元格即可
#############################################################################################################
# 设置全局参数
dir = "sister/result/base" # 根据需要修改，这个是你的sister下的目录，例如base、infercnv等
os.makedirs(dir, exist_ok=True)

############################################################################################################# 
# 读取本研究需要的adata文件
#############################################################################################################
def read_files(filename, pattern):
    from pathlib import Path
    import scanpy as sc

    files = list(Path(filename).rglob(pattern))
    for file in files: print(file)
    adata = sc.concat([sc.read_h5ad(f) for f in files])

    return adata
adata = read_files(filename="/home/students/zhangtao/research/singlecell/database", pattern="*(韩国).h5ad") # 需要根据实际的工作路径和文件名称修改
