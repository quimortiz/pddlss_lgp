
import subprocess
import yaml
import numpy as np
import matplotlib.pyplot as plt

# it is working!!
# experiments should be:
# lgp style
# meta style NO CONFLICTS
# meta style CONFLICTS (but I need this robust!)
# stream style NO REUSE
# stream style REUSE?

# TODO: add more noise to initialization??

# table with obstalce in the middle 

# SAME with 6 objects

# Build a tower with 2 robots -- some reaching, some not.

# problems: 

# 1 robot , 2 objs  -- debug

# 4 PROBLEMS ONLY!!

# 1 robot, 4 objs -- REAL

# 1 robot, 4 objs  with obstacle in the middle

# 2 robots, 6 objs -- tower or all in table?

# 2 robots, 2 objects, require 1 move away (e.g. more logic actions than required)

# MAYBE TWO MORE?

# EXTRA: something with TOOL USE and handover?





# GO with two robots!!




# REAL
algs = ["lgp_style", "stream_style", "meta_style"]
# sizes = [.07, .1, .2, .3]
num_runs = 10

# TEST
# algs = ["lgp_style", "stream_style"]
# algs = ["meta_style"]
# algs = ["stream_style"]
# algs = ["stream_style"]
# sizes = [.3]
# sizes = [.1, .12,  .15, .2, .3]
sizes = [.09 , .1 , .12, .15, .2 , .25, .3, .4] 
# sizes = [.12]
# sizes = [.2, .3]
num_runs = 10

import random
import string
import pathlib

color_map = {"lgp_style" : '#1f77b4',
  "stream_style" : '#ff7f0e',
  "meta" : '#2ca02c'}

  # other defautl colors in matplotlib
  # '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])


def gen_random_id(n: int):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=n))

def get_time_stamp():
    import datetime
    now = datetime.datetime.now()
    return now.strftime("%Y-%m-%d_%H-%M-%S")

def experiment_table_run():
    time_stamp = get_time_stamp()


    summary_all = []
    files_summary = []

    for alg in algs:
        for size in sizes:

            base_path = "results/" + "small_vs_big/" + "table_size_" + str(int( 100*size)) + "/" + alg +  "/" + time_stamp 

            print("A")
            files_local = []
            solved_status = []
            for i in range(num_runs):
                print("B")

                base_path_with_run = base_path + "/r" + str(i).zfill(2)

                pathlib.Path(base_path_with_run).mkdir(parents=True, exist_ok=True)


                cmd = ["./x.exe", "--run_test=small_vs_big_table", "--", "-cfg", "algs/" + alg + ".cfg", "-vis", "0", "-table_size", str(size), "-base_path",  base_path_with_run, "-seed", str(i)] 

                print("running cmd: ", ' '.join(cmd))

                with open(base_path_with_run + "/cmd.txt", "w") as f:
                    f.write(' '.join(cmd))

                id = gen_random_id(6)

                f1 = open(f"/tmp/ss_lgp_stdout_{id}.txt", "w")
                f2 = open(f"/tmp/ss_lgp_stderr_{id}.txt", "w")


                subprocess.run(cmd, stdout=f1, stderr=f2)

                f1.close()
                f2.close()


                D1 = {}
                D2 = {}
                with open( base_path_with_run + "/__benchmark.txt" , "r") as f:
                    D1 = yaml.load(f, Loader=yaml.FullLoader)
                    D1 = D1[1:] # remove first entry
                with open( base_path_with_run + "/__report.txt" , "r") as f:
                    D2 = yaml.load(f, Loader=yaml.FullLoader)

                D2["nlp"] =  D1
                D2["total_time"] = sum([d["ms"] for d in D1])
                D2["num_opts"] = len(D1)

                with open( base_path_with_run + "/__report_all.txt" , "w") as f:
                    yaml.dump(D2, f)

                files_local.append(base_path_with_run + "/__report_all.txt")

                # files_local.append(base_path_with_run + "/__benchmark.txt")
                # solved_status.append(base_path_with_run + "/__report.txt")

            # summary of all runs
            
            Ds = []
            summary = {}
            for file in files_local:
                with open(file, "r") as f:
                    Ds.append(yaml.load(f, Loader=yaml.FullLoader))

            print("Ds: ", Ds)
            average = np.mean([d["total_time"] for d in Ds])
            std_dev = np.std([d["total_time"] for d in Ds])
            print("average: ", average)
            print("std_dev: ", std_dev)

            average = np.mean([d["num_opts"] for d in Ds])
            std_dev = np.std([d["num_opts"] for d in Ds])
            print("average: ", average)
            print("std_dev: ", std_dev)

            summary["time_stamp"] = time_stamp
            summary["success_rate"] = sum(d["solved"] for d in Ds) / len(Ds)
            summary["alg"] = alg
            summary["size"] = size
            summary["files"] = files_local
            summary["time_average"] = float(np.mean([d["total_time"] for d in Ds]))
            summary["time_std"] = float(np.std([d["total_time"] for d in Ds]))
            summary["num_opts_average"] = float(np.mean([d["num_opts"] for d in Ds]))
            summary["num_opts_std"] = float(np.std([d["num_opts"] for d in Ds]))

            summary_file = base_path + "/summary.yaml"
            with open(summary_file, "w") as f:
                yaml.dump(summary, f)

            files_summary.append(summary_file)
            summary_all.append(summary)

    return files_summary

def experiment_table_eval(files):

    print("input files: ", files)
    save_input_file = f"/tmp/experiment_table_input_{gen_random_id(6)}.yaml"
    print("save_input_file: ", save_input_file)
    with open(save_input_file, "w") as f:
        yaml.dump(files, f)


    # print("input files: ", files)
    # with open("/tmp/experiment_table_input.yaml", "r") as f:
    #     files = yaml.safe_load(f)



    Ds =[]

    for file in files:
        with open(file, "r") as f:
            Ds.append(yaml.load(f, Loader=yaml.FullLoader))

    print("Ds: ", Ds)


    success_rate_threshold = .5 
    for alg_name in algs:


        D_alg = [d for d in Ds if d["alg"] == alg_name]


        if len(D_alg) == 0:
            continue

        x = [d["size"] for d in D_alg]
        y = [d["time_average"] if d["success_rate"] > success_rate_threshold else np.inf for d in D_alg]
        error = [d["time_std"] / np.sqrt(len(D_alg)) if d["success_rate"] > success_rate_threshold else np.inf for d in D_alg] 

        # plt.plot(x,y,label=alg_name, marker="o")

        plt.errorbar(x, y, yerr=error, fmt='o', label=alg_name, color=color_map.get(alg_name, "red"))
        # plt.set_title('variable, symmetric error')



    plt.legend()

    plt.ylabel("time [ms]")
    plt.xlabel("table size")

    fileout_fig = f"results_plots/small_vs_big/{get_time_stamp()}.pdf"
    pathlib.Path(fileout_fig).parent.mkdir(parents=True, exist_ok=True)
    with open(fileout_fig  + ".log", "w") as f:
        yaml.dump({"files": files} , f)


    plt.savefig(fileout_fig)
    plt.show()





def experiment_table():
    # SOLVE in LGP STYLE

    files = []
    files = experiment_table_run()

    # with open("/tmp/experiment_table_input.yaml", "r") as f:
    #     files = yaml.safe_load(f)

    experiment_table_eval(files)


experiment_table()


# files = []
# reload = "/tmp/experiment_table_input_X17VIK.yaml"
# with open(reload, "r") as f:
#     files = yaml.safe_load(f)
#
# experiment_table_eval(files)

