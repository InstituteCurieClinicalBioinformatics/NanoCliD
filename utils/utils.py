import subprocess
import os
import glob
import pandas as pd
from snakemake.io import expand
from itertools import product
from nanoclid import NanoClid

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

def getChrom(bed, genomeFile, outFolder, split=False):
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
                df[df[0] == chromosome].to_csv(f"{outFolder}/{os.environ['USER']}_{chromosome}_clair3.bed", sep = "\t", header = False, index = None)
    return chromosomes

def sendMail(config, containersPath, curieNetwork, handler, log=None):
    if handler == "onstart":
        if curieNetwork:
            cmd = f"""export SINGULARITYENV_PDIR={config["workspace"]['pDir']} && singularity exec \
            -B {os.path.dirname(config["workspace"]["pDir"])} \
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
            cmd = f"""export SINGULARITYENV_PDIR={config["workspace"]['pDir']} && singularity exec \
            -B {os.path.dirname(config["workspace"]["pDir"])} \
            --no-home \
            --cleanenv \
            {os.path.join(containersPath, config['bioInfoCliTools']['sif'])} \
            python3 \
            /usr/local/code3/curie/SendMail.py \
            -a 1 \
            -t {config['email']} \
            -s '{config['errorMail']['subject']}' \
            -c '{config['errorMail']['content']}' \
            --custom 1"""
        else:
            cmd = f'mail -s "An error occured for NanoCliD" {config["email"]} < {log}'
    subprocess.call(cmd, shell = True)

#onlyOnePath -> si on est en multiplex, les samples ont tous le meme folder de pod5/fast5, donc on ne renvoie qu'un seul dossier par injection
def getFast5Dir(wildcards, fromBlow5=None, extension=None, onlyOnePath=None, config=None):
    if hasattr(wildcards, "mergedSamples"):
        return config["fast5Dir"][wildcards.mergedSamples]
    if extension and fromBlow5 == "yes":
        return config["fast5Dir"][wildcards.samples].replace('blow5', extension)
    if onlyOnePath:
        #config["samples"] = ['RB01_1', 'RB01_2', 'RB02_1', 'RB02_2', 'RB03_1', 'RB03_2'] par exemple si 2 injections avec 3 samples multiplexes
        barcode = [sample for sample in config["wildcards"]["samples"] if f'_{wildcards.injections}'][0] #RB01_1 ou RB01_2 selon le num d'injection
        return config["fast5Dir"][barcode]
    return config["fast5Dir"][wildcards.samples]

def getReportFile(wildcards, config):
    return config["reportFiles"][wildcards.samples]

def getSubSampling(wildcards, outDir, extension, config):
    path = {"bam" : os.path.join(outDir, "INJECTION", "SAMPLE", "Mapping", "SAMPLE.bam"), \
    "fastq" : os.path.join(outDir, "INJECTION", "SAMPLE", "FASTQ", "SAMPLE.fastq.gz"), \
    "blow5" : os.path.join(outDir, "INJECTION", "SAMPLE", "blow5"), \
    "fast5" : os.path.join(outDir, "INJECTION", "SAMPLE", "fast5"), \
    "modkit" : os.path.join(outDir, 'INJECTION', 'SAMPLE', 'BAM_MODIFIED_MODKIT', 'SAMPLE.bam'), \
    "deepmod" : os.path.join(outDir, 'INJECTION', 'SAMPLE', 'BAM_MODIFIED_DEEPMOD2', 'SAMPLE.bam')}
    files = []
    if hasattr(wildcards, "mergedSamples"):
        samplesToMerge = config["samplesToMerge"][wildcards.mergedSamples]
        for injection in config["wildcards"]["injections"]:
            for sample in config["fast5Dir"].keys():
                if sample in samplesToMerge and injection in config["fast5Dir"][sample]:
                    files.append(path[extension].replace("INJECTION", injection).replace("SAMPLE", sample))
    else:
        files.append(path[extension].replace("INJECTION", wildcards.injections).replace("SAMPLE", wildcards.samples))
    return files

def concatFastq(inputFolder, run, injections, samples, outDir):
    folders = NanoClid.getFiles(os.path.join(inputFolder, run), "d", "fastq_pass")
    if folders != "":
        for folder in folders:
            for sample in samples:
                nb = sample.split("RB")[-1].split("_")[0]
                if os.path.exists(os.path.join(folder, f"barcode{nb}")):
                    for injection in injections:
                        if injection in folder:
                            cmd = f"mkdir -p {outDir}/{injection}/{sample}/FASTQ && \
                            cat {folder}/barcode{nb}/*.fastq.gz > {outDir}/{injection}/{sample}/FASTQ/{sample}.fastq.gz"
                            subprocess.call(cmd, shell = True)


def getClosestBin(cnvFromBamBin):
    binList = [1000, 5000, 10000, 15000, 30000, 50000, 100000, 500000, 1000000, 2500000, 5000000]
    with open(cnvFromBamBin, "r") as cnvFromBamBin:
        cnvFromBamBin = int(cnvFromBamBin.read())
    diffList = [abs(cnvFromBamBin-val) for val in binList]
    return binList[diffList.index(min(diffList))]    
