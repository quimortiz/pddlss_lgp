
import subprocess
import yaml
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List
import pandas
import multiprocessing
from multiprocessing import Pool
import shutil
import os

from matplotlib.backends.backend_pdf import PdfPages

print_lock = multiprocessing.Lock()


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

# 1 robot , 2 objs  -- debug DONE

# 4 PROBLEMS ONLY!!

# 1 robot, 4 objs -- REAL DONE

# 1 robot, 4 objs  with obstacle in the middle

# 2 robots, 6 objs -- tower or all in table?

# 2 robots, 2 objects, require 1 move away (e.g. more logic actions than
# required)

# MAYBE TWO MORE?

# EXTRA: something with TOOL USE and handover?


# GO with two robots!!


n_cores = -1


if n_cores == -1:
    n_cores = int(multiprocessing.cpu_count() / 2)

# REAL


@dataclass
class Problem:
    configC: str
    goal: str
    configL: str
    name_base: str
    opts: List[str] = field(default_factory=list)

    def args(self):
        return [
            "-confFile",
            f"\"{self.configC}\"",
            "-goalFile",
            f"\"{self.goal}\"",
            "-folFile",
            f"\"{self.configL}\""]

    def additonal_args(self):
        return self.opts

    def name(self):
        return self.configC.split("/")[-1].split(".")[0]


problems = [

    Problem(configC="./models/table_experiment.g",
            goal="./models/table_experiment_goal.g",
            configL="./fol_lab_bench_auto.g",
            name_base="small_vs_big2",
            opts=["-ct/max_depth_discrete", "4"]),


    Problem(configC="./models/table_experiment_4.g",
            goal="./models/table_experiment4_goal.g",
            configL="./fol_lab_bench_auto.g",
            name_base="small_vs_big2",
            opts=["-ct/max_depth_discrete", "8"]),

    Problem(configC="./models/table_experiment_4_obstacle.g",
            goal="./models/table_experiment4_goal.g",
            configL="./fol_lab_bench_auto.g",
            name_base="small_vs_big4_obs",
            opts=["-ct/max_depth_discrete", "8"])

]


problems_experiment2 = [
    "experiment2/problem1.cfg",
    "experiment2/problem2.cfg",
    "experiment2/problem3.cfg",
    # "experiment2/problem4.cfg"
]

problems_experiment2_easy = [
    "experiment3/problem2.cfg",
    "experiment3/problem3.cfg",
    "experiment3/problem4.cfg",
    "experiment3/problem5.cfg",
    # "experiment2/problem3.cfg",
    # "experiment2/problem4.cfg"
]


algs_experiment2 = [
    "experiment2/algs/meta_w_conflicts.cfg",
    "experiment2/algs/meta_wo_conflicts.cfg",
    "experiment2/algs/stream.cfg",
    "experiment2/algs/lgp.cfg"
]


problems = [problems[-1]]


# algs = ["stream_style"]
# algs = ["lgp_style", "stream_style", "meta_style"]
algs = ["lgp_style"]
# sizes = [.07, .1, .2, .3]

# TEST
# algs = ["lgp_style", "stream_style"]
# algs = ["meta_style"]
# algs = ["stream_style"]
# algs = ["stream_style"]
# sizes = [.3]
# sizes = [.1, .12,  .15, .2, .3]
# sizes = [.09 , .1 , .12, .15, .2 , .25, .3, .4]
# sizes = [.12]
# sizes = [.2, .3]

# sizes = [.3]
# sizes = [.2,.25,.3,.4,.5]
sizes = [.4]
# sizes = [.15]
num_runs = 1

import random
import string
import pathlib
from pathlib import Path

color_map = {"lgp_style": '#1f77b4',
             "stream_style": '#ff7f0e',
             "meta": '#2ca02c',
             "xx": '#d62728'}

# other defautl colors in matplotlib
# '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])


def gen_random_id(n: int = 6):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=n))


def get_time_stamp():
    import datetime
    now = datetime.datetime.now()
    return now.strftime("%Y-%m-%d_%H-%M-%S")


def run_cmd(cmd: str) -> None:

    id = gen_random_id()
    filename_stdout = f"/tmp/ss_lgp_stdout_{id}.txt"
    fileanem_stderr = f"/tmp/ss_lgp_stderr_{id}.txt"

    print_lock.acquire()
    print("running cmd: ", ' '.join(cmd))
    print("filename_stdout: ", filename_stdout)
    print("fileanem_stderr: ", fileanem_stderr)
    print_lock.release()

    f1 = open(filename_stdout, "w")
    f2 = open(fileanem_stderr, "w")

    subprocess.run(cmd, stdout=f1, stderr=f2)

    f1.close()
    f2.close()

    D1 = {}
    D2 = {}
    base_path_with_run = cmd[-3]

    __file = base_path_with_run + "/__benchmark.txt"

    if not pathlib.Path(__file).is_file():
        raise Exception("file not found: ", __file)

    with open(__file, "r") as f:
        D1 = yaml.load(f, Loader=yaml.FullLoader)
        D1 = D1[1:]  # remove first entry

    __file = base_path_with_run + "/__report.txt"

    # RUN CMND
    if not pathlib.Path(__file).is_file():
        raise Exception("file not found: ", __file)

    with open(base_path_with_run + "/__report.txt", "r") as f:
        D2 = yaml.load(f, Loader=yaml.FullLoader)

    D2["nlp"] = D1
    D2["total_time"] = sum([d["ms"] for d in D1])
    D2["num_opts"] = len(D1)

    with open(base_path_with_run + "/__report_all.txt", "w") as f:
        yaml.dump(D2, f)


def experiment_run(problems: List[str], algs: List[str], num_runs: int = 1):
    time_stamp = get_time_stamp()

    summary_all = []
    files_summary = []
    cmds = []
    experiments = []
    index = 0

    for problem, alg in [(p, a) for p in problems for a in algs]:

        alg_short = alg.split("/")[-1].split(".")[0]

        base_path = "results/" + \
            problem.split(".")[0] + "/" + alg_short + "/" + time_stamp
        print("base_path: ", base_path)

        files_local = []
        for i in range(num_runs):

            base_path_with_run = base_path + "/r" + str(i).zfill(2) + "/"

            pathlib.Path(base_path_with_run).mkdir(parents=True, exist_ok=True)

            cmd = [
                "./x.exe",
                "--run_test=small_vs_big_table",
                "--",
                "-cfg",
                problem,
                "-cfg2",
                alg,
                "-base_path",
                base_path_with_run,
                "-seed",
                str(i)]

            with open(base_path_with_run + "/cmd.txt", "w") as f:
                f.write(' '.join(cmd))
                f.write('\n')

            cmds.append(cmd)
            files_local.append(base_path_with_run)
        experiments.append((problem, alg, files_local, base_path))

    print(f"Start a pool with {n_cores}:")

    with Pool(n_cores) as p:
        p.map(run_cmd, cmds)

        # summary of all runs

    # TODO: I could even parallelize this!
    for data in experiments:
        problem, alg, files_local, base_path = data
        Ds = []
        summary = {}
        for file in files_local:
            with open(file + "__report_all.txt", "r") as f:
                Ds.append(yaml.load(f, Loader=yaml.FullLoader))

        print("Ds: ", Ds)
        average = np.median([d["total_time"] for d in Ds if d["solved"]])
        std_dev = np.std([d["total_time"] for d in Ds if d["solved"]])
        print("average: ", average)
        print("std_dev: ", std_dev)

        average = np.median([d["num_opts"] for d in Ds if d["solved"]])
        std_dev = np.std([d["num_opts"] for d in Ds if d["solved"]])
        print("average: ", average)
        print("std_dev: ", std_dev)

        percentil_up = .9
        percentil_down = .1

        out = np.percentile([d["total_time"] for d in Ds], [
                            percentil_down * 100, percentil_up * 100])

        # np.percentile(
        summary["time_stamp"] = time_stamp
        summary["success_rate"] = sum(d["solved"] for d in Ds) / len(Ds)
        summary["alg"] = alg
        summary["problem"] = problem
        summary["files"] = files_local
        summary["time_percentil_down"] = float(out[0])
        summary["time_percentil_up"] = float(out[1])
        summary["time_average"] = float(
            np.median([d["total_time"] for d in Ds]))
        summary["time_std"] = float(np.std([d["total_time"] for d in Ds]))
        summary["num_opts_average"] = float(
            np.mean([d["num_opts"] for d in Ds]))
        summary["num_opts_std"] = float(np.std([d["num_opts"] for d in Ds]))

        summary_file = base_path + "/summary.yaml"
        with open(summary_file, "w") as f:
            yaml.dump(summary, f)

        files_summary.append(summary_file)
        summary_all.append(summary)

    # store the results as the csv
    # TODO
    print("summary_all: ", summary_all)

    df = pandas.DataFrame(summary_all)
    fileout = f"results/summary_{time_stamp}.csv"
    print("writing to file: ", fileout)
    df.to_csv(fileout)

    fileout = f"results/summary_{time_stamp}.yaml"
    print("writing to file: ", fileout)
    with open(fileout, "w") as f:
        yaml.dump(summary_all, f)

    fileout = f"results/summary_{time_stamp}.yaml.log"
    print("writing to file: ", fileout)
    with open(fileout, "w") as f:
        yaml.dump(files_summary, f)

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

    Ds = []

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
        y = [d["time_average"] if d["success_rate"] >
             success_rate_threshold else np.inf for d in D_alg]
        error = [
            d["time_std"] /
            np.sqrt(
                len(D_alg)) if d["success_rate"] > success_rate_threshold else np.inf for d in D_alg]

        # plt.plot(x,y,label=alg_name, marker="o")

        plt.errorbar(
            x,
            y,
            yerr=error,
            fmt='o',
            label=alg_name,
            color=color_map.get(
                alg_name,
                "red"))
        # plt.set_title('variable, symmetric error')

    plt.legend()

    plt.ylabel("time [ms]")
    plt.xlabel("table size")

    fileout_fig = f"results_plots/small_vs_big/{get_time_stamp()}.pdf"
    pathlib.Path(fileout_fig).parent.mkdir(parents=True, exist_ok=True)
    with open(fileout_fig + ".log", "w") as f:
        yaml.dump({"files": files}, f)

    plt.savefig(fileout_fig)
    # plt.show()


def tests():
    # TODO
    cmd = "m4 && gdb --args  ./x.exe --run_test=small_vs_big_table -- -cfg test_plans/test_two_robots.cfg -vis | tee quim.txt"


def generate_texpdf(filename_tex: str) -> None:
    lines = [
        r"\documentclass{standalone}",
        r"\usepackage{amsmath}",
        r"\usepackage{amsfonts}",
        r"\usepackage{siunitx}",
        r"\usepackage{booktabs}", 
        r"\usepackage{xcolor}", 
        r"\begin{document}",
        r"\input{" + filename_tex + "}",
        r"\end{document}\n",
    ]

    print("writing to ", filename_tex + ".tex")
    with open(filename_tex + ".tex", "w") as f:
        f.write('\n'.join(lines))

    pathlib.Path("/tmp/dynoplan/").mkdir(parents=True, exist_ok=True)
    pathlib.Path("/tmp/dynoplan/tex").mkdir(parents=True, exist_ok=True)
    f_stdout = open("/tmp/dynoplan/latexmk_out.log", "w")
    f_stderr = open("/tmp/dynoplan/latexmk_err.log", "w")

    out_dir = "/tmp/dynoplan/tex/"
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    cmd = f"latexmk -f -pdf -output-directory={out_dir} -interaction=nonstopmode {filename_tex}.tex".split()

    print("running cmd: ", " ".join(cmd))
    subprocess.run(cmd, stdout=f_stdout, stderr=f_stderr)

    f_stdout.close()
    f_stderr.close()

    print("latexmk stdout:")
    os.system("cat /tmp/dynoplan/latexmk_out.log")
    print("latexmk stderr:")
    os.system("cat /tmp/dynoplan/latexmk_err.log")

    p = Path(filename_tex)
    pdf_name = "/tmp/dynoplan/tex/" + str(p.stem) + ".tex.pdf"

    print(f"copy  {pdf_name} to {filename_tex + '.pdf'}")
    shutil.copy(pdf_name, filename_tex + ".pdf")

    tmp_location = "/tmp/dynoplan/table_tex.pdf"

    print(f"copy  {pdf_name} to {tmp_location}")
    shutil.copy(pdf_name, tmp_location)



def experiments_eval(files):
    """
    """
    print("input files", files)

    fileout_input = f"/tmp/ss_{gen_random_id()}.csv"
    print(f"fileout_input {fileout_input}")
    with open(fileout_input, "w") as f:
        yaml.dump(files, f)

    data = []
    for file in files:

        with open(file, "r") as f:
            data.append(yaml.safe_load(f))

    print("data is ")
    print(data)

    # FOR NOW: I only have one problem!!
    success_rate_threshold = .5

    # how many problems?

    problems = list(set([d["problem"] for d in data]))
    algs = list(set([d["alg"] for d in data]))

    color_map = {"lgp_style": '#1f77b4',
                 "stream_style": '#ff7f0e',
                 "meta": '#2ca02c',
                 "xx": '#d62728'}

    def alg_to_color(alg: str) -> str:
        out = ""
        if "lgp" in alg:
            out = "red"
        elif "stream" in alg:
            out = "blue"
        elif "meta_wo_conflicts" in alg:
            out = "green"
        elif "meta_w_conflicts" in alg:
            out = "orange"
        else:
            out = "black"
        return out

    problems.sort()
    print("problems", problems)
    # sys.exit(0)

    x = [i for i in range(len(problems))]

    fileout_fig = f"results_plots/transfer/{get_time_stamp()}.pdf"
    pathlib.Path(fileout_fig).parent.mkdir(parents=True, exist_ok=True)

    pp = PdfPages(fileout_fig)

    # get the best for each problem
    

    problem2bestTime = {}
    for problem in problems:
            
        D_prob = [d for d in data if d["problem"] == problem]

        # best alg for this problem
        D_prob.sort(key=lambda x: x["time_average"])

        best_time = D_prob[0]["time_average"]
        problem2bestTime[problem] = best_time


    print("problem2bestTime")
    print(problem2bestTime)


    # create a table

    """
                                    ALG 1           |      ALG 2     |     ALG 3 
                              worse | mean | best                     

    problemID | problem | 
    """


    def escape_underscore(s: str) -> str:
        return s.replace("_", "\\_")



    num_data_cols = 3 * len(algs)


    header_line =  f"\\begin{{tabular}}{{||l|l|" + num_data_cols * (  "|" + 3 * "c|" ) + "|}\n"



    algs_names = [ alg.split("/")[-1].split(".")[0] for alg in algs ]

    line1a = ["\\#", "Problem" ] +  [ f"\\multicolumn{{3}}{{c||}}{{{escape_underscore(alg)}}}" for alg in algs_names ] 
    line1b = ["", ""] + len(algs) * ( [ "best" , "median"  , "worse" ])


    # algs
    # line1 = ["problemID", "problem" ] + [ escape_underscore(i+j) for i in algs for j in ["worse", "mean", "best"] ]
    lines_data = []

    def is_highest(val, field: str, algs) -> bool:
        for i, a in enumerate(algs):
            vv = [ d for d in data if d["alg"] == a and d["problem"] == p]
            assert(vv)
            v = vv[0]
            if val < v[field] :
                return False
        return True

    def is_lowest(val, field: str, algs) -> bool:
        for i, a in enumerate(algs):
            vv = [ d for d in data if d["alg"] == a and d["problem"] == p]
            assert(vv)
            v = vv[0]
            if val > v[field] :
                return False
        return True


    for i, p in enumerate(problems):
        line = []
        line += [i, escape_underscore(p)]
        for a  in algs:
            vv = [ d for d in data if d["alg"] == a and d["problem"] == p]

            assert(vv)
            v = vv[0]
            for f in ["time_percentil_down", "time_average", "time_percentil_up"]:
                # check if the time_percentil is down

                # CONTINUE HERE with red on the worse 
                # And Green on the the best median.
                # also add succes rate.
                # is_highest = is_worse(f)
                if f == "time_percentil_up" and  is_highest(v[f], f , [ b for b in algs if b != a ] ):
                    line.append(f"\\textcolor{{red}}{{{v[f]:.0f}}}")
                elif f == "time_percentil_up" and is_lowest(v[f], f , [ b for b in algs if b != a ] ):
                    line.append(f"\\textcolor{{green}}{{{v[f]:.0f}}}")
                elif f == "time_average" and is_lowest(v[f], f , [ b for b in algs if b != a ] ):
                    line.append(f"\\textcolor{{green}}{{{v[f]:.0f}}}")
                else:
                    line.append(f"{v[f]:.0f}")
        lines_data.append(line)

    # lets print some latex!
    

    filename_tex = "/tmp/table.tex"

    with open(filename_tex, "w") as f:
        # f.write(f"\\begin{{tabular}}{{|l|l|{num_data_cols * 'c|'}}}\n")
        f.write(header_line)
        # f.write(" & ".join(line1) + "\\\\\n")
        f.write(" & ".join(line1a) + "\\\\\n")

        f.write(f"\\hline\n")
        f.write(" & ".join(line1b) + "\\\\\n")
        f.write(f"\\hline\n")
        for line in lines_data:
            f.write( " & ".join([str(x) for x in line]) + "\\\\\n")
            f.write(f"\\hline\n")
        f.write(f"\\end{{tabular}}\n")

    generate_texpdf(filename_tex)




    for aa, alg_name in enumerate(algs):

        D_alg = [d for d in data if d["alg"] == alg_name]
        y = []
        error = []
        percentil_ups = []
        percentil_downs = []
        for p in problems:
            D_algprob = [d for d in D_alg if d["problem"] == p]
            assert len(D_algprob) == 1
            d = D_algprob[0]
            xi = d["time_average"] if d["success_rate"] > success_rate_threshold else np.inf
            y.append(xi / problem2bestTime[p])
        # y = [d["time_average"] if d["success_rate"] >
        # success_rate_threshold else np.inf for p in problems in D_algprob]
            ei = d["time_std"] if d["time_std"] > success_rate_threshold else np.inf
            ei = ei / np.sqrt(len(D_algprob))
            error.append(ei)
            percentil_up = d["time_percentil_up"]
            percentil_down = d["time_percentil_down"]

            percentil_ups.append(percentil_up / problem2bestTime[p])
            percentil_downs.append(percentil_down / problem2bestTime[p])

        # plt.errorbar(
        #     x,
        #     y,
        #     yerr=error,
        #     fmt='o',
        #     capsize=10,
        #     label=alg_name,
        #     color=alg_to_color(alg_name))

        plt.errorbar(
            x + aa * .05 * np.ones(len(x)),
            y,
            # yerr=error,
            yerr=[percentil_downs, percentil_ups],
            fmt='o',
            capsize=10,
            label=alg_name)

    plt.legend()
    plt.ylabel("time [ms]")
    plt.xlabel("problem ")
    plt.xticks(x, [problem.split('/')[1].split('.')[0]
               for problem in problems])

    plt.xticks(rotation='vertical')

    # , rotation='vertical

    with open(fileout_fig + ".log", "w") as f:
        yaml.dump({"files": files}, f)

    print("save to file: ", fileout_fig)
    print("save to file: ", fileout_fig + ".log")

    plt.tight_layout()

    plt.savefig(fileout_fig)

    tmp_copy = "/tmp/pddlss_gnlp.pdf"
    print("save to file: ", tmp_copy)
    shutil.copyfile(fileout_fig, tmp_copy)

    plt.clf()

    # ONLY Succes rate
    for aa, alg_name in enumerate(algs):

        D_alg = [d for d in data if d["alg"] == alg_name]
        y = []
        error = []
        percentil_ups = []
        percentil_downs = []
        for p in problems:
            D_algprob = [d for d in D_alg if d["problem"] == p]
            assert len(D_algprob) == 1
            d = D_algprob[0]
            xi = d["success_rate"]
            y.append(xi)
        # y = [d["time_average"] if d["success_rate"] >
        # success_rate_threshold else np.inf for p in problems in D_algprob]
            # ei = d["time_std"] if d["time_std"] > success_rate_threshold else np.inf
            # ei = ei / np.sqrt(len(D_algprob))
            # error.append(ei)
            # percentil_up = d["time_percentil_up"]
            # percentil_down = d["time_percentil_down"]
            #
            # percentil_ups.append(percentil_up)
            # percentil_downs.append(percentil_down)

        # plt.errorbar(
        #     x,
        #     y,
        #     yerr=error,
        #     fmt='o',
        #     capsize=10,
        #     label=alg_name,
        #     color=alg_to_color(alg_name))

        plt.plot(x + aa * .05 * np.ones(len(x)), y, 'o', label=alg_name)

    plt.legend()
    fileout2 = fileout_fig + ".rate.pdf"
    print("save to file: ", fileout2)
    plt.savefig(fileout2)

    tmp_copy = "/tmp/pddlss_gnlp_rate.pdf"
    print("save to file: ", tmp_copy)
    shutil.copyfile(fileout2, tmp_copy)

    # plt.show()


def experiment_transfer():
    # SOLVE in LGP STYLE

    # problems = problems_experiment2
    # algs = algs_experiment2
    # files = []
    num_runs = 3

    problems = problems_experiment2_easy
    algs = algs_experiment2
    files = experiment_run(problems, algs, num_runs)

    # filename =  "results/summary_2023-08-28_09-16-15.yaml.log"

    # filename = "results/summary_2023-08-28_09-20-48.yaml.log"
    # 01-joint-vs-conditional/results/summary_2023-08-27_16-55-04.yaml.log

    #

    # filename = "results/summary_2023-08-28_11-27-11.yaml.log"
    #
    # with open(filename, "r") as f:
    #     files = yaml.safe_load(f)

    experiments_eval(files)


problems_experiment5 = ["experiment5/problem1.cfg",
                        "experiment5/problem2.cfg"]


def experiment_move_block_away():
    num_runs = 10
    problems = problems_experiment5
    algs = algs_experiment2

    algs = ["experiment2/algs/meta_w_conflicts.cfg",
            # "experiment2/algs/meta_wo_conflicts.cfg",
            # "experiment2/algs/stream.cfg"
            ]

    files = experiment_run(problems, algs, num_runs)
    experiments_eval(files)


def experiment_table():
    # SOLVE in LGP STYLE

    files = []

    number_blocks: int = 2
    middle_obstacle: bool = True
    problems: List[str] = []
    stack: bool = True
    if number_blocks == 2:
        if middle_obstacle:
            if stack:
                problems = [
                    "experiment1/problem1o_stack_10.cfg",  # LGP is worse here
                    "experiment1/problem1o_stack_20.cfg"]

            else:
                problems = [
                    "experiment1/problem1o_20.cfg",  # I can use this one: lgp is worse
                    "experiment1/problem1o_30.cfg"
                ]

        else:
            problems = [
                "experiment1/problem1_09.cfg",
                "experiment1/problem1_10.cfg",  # I can use this one: meta is worse
                "experiment1/problem1_15.cfg",
                "experiment1/problem1_20.cfg"
            ]

    if number_blocks == 4:
        if middle_obstacle:
            problems = [
                "experiment1/problem3_25.cfg",
                "experiment1/problem3_30.cfg",
                "experiment1/problem3_40.cfg"
            ]

        else:
            problems = [
                # "experiment1/problem2_15.cfg",
                "experiment1/problem2_20.cfg",
                "experiment1/problem2_30.cfg"
            ]
        # continue here!!

    algs = [
        "algs/lgp_style.cfg",
        "algs/stream_style.cfg",
        # "algs/meta_style.cfg",
        # "algs/meta_style_old.cfg",
        # "algs/meta_style_rand.cfg"
        # meta_wo_conflicts.cfg"
    ]

    num_runs = 5

    files = experiment_run(problems, algs, num_runs)

    # file = "/tmp/ss_H0FUC8.csv"
    #
    #
    # with open(file, "r") as f:
    #     files = yaml.safe_load(f)

    experiments_eval(files)

    # k
    # with open("/tmp/experiment_table_input.yaml", "r") as f:
    #     files = yaml.safe_load(f)

    # experiment_table_eval(files)


# experiment_table()


def experiment_table_more():

    algs = [
        "algs/lgp_style.cfg",
        "algs/stream_style.cfg",
        # "algs/meta_style.cfg",
        # "algs/meta_style_old.cfg",
        # "algs/meta_style_rand.cfg"
        # meta_wo_conflicts.cfg"
    ]

    num_runs = 10
    problems = [
        # "experiment1/transfer_try1.cfg",
        # "experiment1/handover_try1.cfg",
        "experiment1/handover_try1.cfg",
        "experiment1/handover_try1_25.cfg",  # OPTIMIZATION is much better here!!
    ]

    # files = "/tmp/ss_G6O5PI.csv"
    # #
    # with open(files, "r") as f:
    #     files = yaml.safe_load(f)

    files = experiment_run(problems, algs, num_runs)

    experiments_eval(files)


def experiment_move_away():
    # Simple experiment to move away blocks,
    # using one robot.

    # What should I use?

    algs = [
        "algs/lgp_style.cfg",
        "algs/stream_style.cfg",
        "algs/meta_style.cfg",
        "algs/meta_style_old.cfg",
        "algs/meta_style_rand.cfg"
        # meta_wo_conflicts.cfg"
    ]

    num_runs = 10
    problems = [
        "experiment1/problem_move_away_1_small.cfg",
        # LGP is better :) -- seems that we have results :)
        "experiment1/problem_move_away_1.cfg"
    ]

    # files = "/tmp/ss_G6O5PI.csv"
    # #
    # with open(files, "r") as f:
    #     files = yaml.safe_load(f)

    files = experiment_run(problems, algs, num_runs)

    experiments_eval(files)


# experiment_table_more()

# experiment_move_away()

problems_thesis = [
    "experiment1/problem_move_away_1_small.cfg",
    "experiment1/handover_try1_25.cfg", 
    "experiment1/transfer_try1.cfg",
    "experiment1/problem1o_stack_10.cfg",
    "experiment1/problem1o_20.cfg", 
    "experiment1/problem1_10.cfg",  
]

algs_thesis = [
        "algs/lgp_style.cfg",
        "algs/stream_style.cfg",
        "algs/meta_style.cfg",
        "algs/meta_style_old.cfg",
        "algs/meta_style_rand.cfg"
]


def experiments_thesis():
    # num_runs = 10
    # files = experiment_run(problems_thesis, algs_thesis, num_runs)

    files_log = "results_plots/transfer/2023-09-06_11-27-29.pdf.log"                                                                                                                                                      

    with open(files_log) as f:
        files = yaml.safe_load(f)["files"]

    experiments_eval(files)

experiments_thesis()




# experiment_transfer()

# experiment_move_block_away()


# experiments_two_robots_three_tables()


# files = []
# reload = "/tmp/experiment_table_input_X17VIK.yaml"
# with open(reload, "r") as f:
#     files = yaml.safe_load(f)
#
# experiment_table_eval(files)
