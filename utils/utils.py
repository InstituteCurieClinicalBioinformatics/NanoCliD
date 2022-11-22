import subprocess
import os
from snakemake.io import expand

def getFiles(path, kind, pattern, files=False):
    find = subprocess.check_output(f"find {path} -type {kind} -name {pattern}").decode("utf-8").rstrip().split("\n")
    if len(find) == 0:
        print(f"No files found in path {path} with pattern {pattern}")
        return ""
    elif files:
        return find
    else:
        return find[0]

def getSamplesDir(run):
    folders = os.listdir(run)
    samples = []
    for folder in folders:
        if os.path.isdir(folder):
            samples.append(folder)
    return samples

def getOutputs(template, RES_FOLDER, samples, run_id):
    outputs = []
    with open(template, "r") as f:
        for line in f:
            if "{injection}" in line.rstrip() or "{run}" in line.rstrip():
                outputs.append(expand(os.path.join(RES_FOLDER, run_id, line.rstrip()), injection = samples, run = run_id))
            else:
                outputs.append(os.path.join(RES_FOLDER, run_id, line.rstrip()))
    return outputs