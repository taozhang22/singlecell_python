def create_code_file(wd, dir):
    import os
    os.chdir(wd)

    os.makedirs(f"{dir}/raw")
    os.makedirs(f"{dir}/log")
    os.makedirs(f"{dir}/output")
create_database_file(wd="/home/students/zhangtao/research/singlecell", dir="database")
