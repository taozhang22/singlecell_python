conda activate scenv
nohup papermill "GSE200997(印第安纳大学)/GSE200997(印第安纳大学).ipynb" \
  "output/GSE200997(印第安纳大学)_output.ipynb" \
  -k scenv \
  > "log/GSE200997(印第安纳大学).log" 2>&1 &
