############################################################################################################# 
# 创建scenv环境
############################################################################################################# 
conda activate base
conda create -n scenv python=3.10.14 -y
conda activate scenv
conda install scanpy -y
pip install ipykernel
python -m ipykernel install --user --name=scenv --display-name "scenv"
pip install papermill
pip install ipython-autotime
pip install openpyxl
pip install igraph
pip install leidenalg
pip install scikit-image
pip install scvi-tools[cuda]
pip install harmonypy
pip install infercnvpy
pip install gseapy
pip install sccoda
pip install liana
pip install cellphonedb
