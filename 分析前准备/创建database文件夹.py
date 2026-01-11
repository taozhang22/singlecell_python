def create_database_file(wd, dir):
    import os
    os.chdir(wd)

    os.makedirs(f"{dir}/raw", exist_ok=True)
    os.makedirs(f"{dir}/log", exist_ok=True)
    os.makedirs(f"{dir}/output", exist_ok=True)
create_database_file(wd="/home/students/zhangtao/research/singlecell", dir="database")
