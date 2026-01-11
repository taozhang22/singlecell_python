# environmet: scenv
def prepare_environment(dir):
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
    import scvi
    import singlecell as scell
    os.chdir("/home/students/zhangtao/research/singlecell/python/Sister")

    print(f"Current python version: {sys.version}")
    !conda list
    pd.set_option("display.width", 1000)
prepare_environment(dir="/home/students/zhangtao/research/singlecell/python/Sister")
