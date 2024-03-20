import subprocess
import os
import glob
import pandas as pd
from snakemake.io import expand
from itertools import product

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

def filter_combinator(combinator, keepCombinaisons):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            # Use frozenset instead of tuple
            # in order to accomodate
            # unpredictable wildcard order
            if frozenset(wc_comb) in keepCombinaisons:
                yield wc_comb
    return filtered_combinator

def getOutputs(template, RES_FOLDER, RES_FOLDER_INJECTION, wildcards, combinaison):
    if combinaison != "":
        keepCombinaisons = {frozenset(wc_comb.items()) for wc_comb in combinaison}
        filtered_product = filter_combinator(product, keepCombinaisons)
    wildcardsList = list(wildcards.keys())
    outputs = []
    with open(template, "r") as f:
        for line in f:
            if "{injections}" in line or "{samples}" in line:
                outputs.append(expand(os.path.join(RES_FOLDER_INJECTION, line.rstrip()), filtered_product, **wildcards))
            elif any(wildcard in line for wildcard in wildcardsList):
                outputs.append(expand(os.path.join(RES_FOLDER, line.rstrip()), **wildcards))
            else:
                outputs.append(os.path.join(RES_FOLDER, line.rstrip()))
    return outputs

def get_mem_depending_of_input(wildcards, input, attempt, min_ram_gb=None):
    return min_ram_gb + ((input.size//10000) * attempt * 1.8)

def splitBedForMethylation(bed, genomeFile, chrom=False):
    if os.path.basename(bed) == "empty.bed":
        bed = genomeFile
        chrom = True
    bed = pd.read_csv(bed, sep = "\t", header = None)
    if chrom:
        return list(set(bed[0]))
    else:
        bed[1] = [f'{value:,}' for value in bed[1].tolist()]
        bed[2] = [f'{value:,}' for value in bed[2].tolist()]
        bed[4] = bed[0] + ":" + bed[1] + "-" + bed[2]
        return bed[4].tolist()

def getChrom(bed, genomeFile, split=False):
    if os.path.basename(bed) == "empty.bed":
        bed = genomeFile
        df = pd.read_csv(bed, sep = "\t", header = None)
        df[2] = [0] * df.shape[0]
        df = df[[0,2,1]]
    else:
        df = pd.read_csv(bed, sep = "\t", header = None)
    chromosomes = list(set(df[0]))
    if split:
        for chromosome in chromosomes:
                df[df[0] == chromosome].to_csv(f"/tmp/{chromosome}_clair3.bed", sep = "\t", header = False, index = None)
    return chromosomes

    shell(f"""singularity  exec --no-home --cleanenv {os.path.join(CONTAINERS_PATH, config['bioInfoCliTools']['sif'])} \
            python3 /usr/local/code3/curie/SendMail.py -a 1 -t {config['errorMail']['bic']} \
            -s "[NanoCliD] Run {config['run']} started" \
            -c "" --custom 1""")   

def sendMail(config, containersPath, curieNetwork, handler, log=None):
    if handler == "onstart":
        if curieNetwork:
            cmd = f"""singularity exec \
            --no-home \
            --cleanenv \
            {os.path.join(containersPath, config['bioInfoCliTools']['sif'])} \
            python3 \
            /usr/local/code3/curie/SendMail.py \
            -a 1 \
            -t {config['email']} \
            -s "[NanoCliD] Run {config['run']} started" \
            -c "" \
            --custom 1"""
        else:
            cmd = f"""echo | mail -s "[NanoCliD] Run {config['run']} started" {config['email']}"""
    else:
        if curieNetwork:
            cmd = f"""singularity exec \
            --no-home \
            --cleanenv \
            {os.path.join(containersPath, config['bioInfoCliTools']['sif'])} \
            python3 \
            /usr/local/code3/curie/SendMail.py \
            -a 1 \
            -t {config['email']} \
            -s {config['errorMail']['subject']} \
            -c {config['errorMail']['content']} \
            --custom 1"""
        else:
            cmd = f'mail -s "An error occured for NanoCliD" {config["email"]} < {log}'
    subprocess.call(cmd, shell = True)