############################################################################################################# 
# 创建代码文件夹
#############################################################################################################
# 参数，根据实际情况修改
wd = "/home/yang/research/bioinformation/singlecell/python"
dir = "practice"

import os
os.chdir(wd)
os.makedirs(f"{dir}/Script/raw", exist_ok=True)
os.makedirs(f"{dir}/Script/log", exist_ok=True)
os.makedirs(f"{dir}/Script/output", exist_ok=True)

############################################################################################################# 
# 代码开头
#############################################################################################################
# environmet: scenv
import time
start = time.time()

import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
os.chdir("/home/students/zhangtao/research/singlecell/python/sister") # 需要根据实际的工作路径修改

print(f"Current python version: {sys.version}")
!conda list
pd.set_option("display.width", 1000)

dir = "result/base" # 根据需要修改，这个是你的sister下的目录，例如base、infercnv等
os.makedirs(dir, exist_ok=True)

############################################################################################################# 
# 代码结尾
#############################################################################################################
end = time.time()
time_elapsed = end - start
print(f"Time elapsed: {time_elapsed/60:.2f} minutes")
