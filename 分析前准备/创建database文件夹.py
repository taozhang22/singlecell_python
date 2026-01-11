def create_code_file(wd, dir):
    import os
    os.chdir(wd)

    os.makedirs("{dir}/raw")
    os.makedirs("{dir}/log")
    os.makedirs("{dir}/output")
create_code_file(wd="/home/students/zhangtao/research/singlecell", dir="database")
