# environmet: scenv
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
