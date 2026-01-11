def read_files(filename, pattern):
    from pathlib import Path
    import scanpy as sc

    files = list(Path(filename).rglob(pattern))
    for file in files: print(file)
    adata = sc.concat([sc.read_h5ad(f) for f in files])

    return adata
adata = read_files(filename="/home/students/zhangtao/research/singlecell/database", pattern="*(韩国)_2.h5ad$")
