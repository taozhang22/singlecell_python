def create_code_file(wd, dir):
    import os
    os.chdir(wd)

    os.makedirs(f"{dir}/Script/raw")
    os.makedirs(f"{dir}/Script/log")
    os.makedirs(f"{dir}/Script/output")
create_code_file(wd="/home/students/zhangtao/research/singlecell/python", dir="Sister")
