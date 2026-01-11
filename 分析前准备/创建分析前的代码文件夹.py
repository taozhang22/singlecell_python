def create_code_file(wd, dir):
    import os
    os.chdir(wd)

    os.makedirs(f"{dir}/Script/raw", exist_ok=True)
    os.makedirs(f"{dir}/Script/log", exist_ok=True)
    os.makedirs(f"{dir}/Script/output", exist_ok=True)
create_code_file(wd="/home/students/zhangtao/research/singlecell/python", dir="sister")
