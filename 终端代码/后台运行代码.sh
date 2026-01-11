conda activate scenv
nohup papermill "/home/students/zhangtao/research/singlecell/python/sister/Script/raw/base.ipynb" \
  "/home/students/zhangtao/research/singlecell/python/sister/Script/output/base_output.ipynb" \
  -k scenv \
  > "/home/students/zhangtao/research/singlecell/python/sister/Script/log/base.log" 2>&1 &
