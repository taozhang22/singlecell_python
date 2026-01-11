def create_code_file(wd, dir):
    import os
    os.chdir(wd)

    os.makedirs("{dir}/Script/raw")
    os.makedirs("{dir}/Script/log")
    os.makedirs("{dir}/Script/output")
create_code_file(wd="/home/students/zhangtao/research/singlecell/python", filename="Sister")
