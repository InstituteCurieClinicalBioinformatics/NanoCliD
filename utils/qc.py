import subprocess
import pandas as pd
from argparse import ArgumentParser
import os
import re
import math

def parseBamStats(bamStatFile):
    stats = []
    keywords = ["average length"]
    for keyword in keywords:
        stats.append(subprocess.check_output(f"grep '{keyword}' {bamStatFile}", shell=True).decode("utf-8").rstrip().split("\t")[-1])
    return stats

def parseMosDepthSummary(mosdepthFile, bedFile):
    df = pd.read_csv(mosdepthFile, sep="\t")
    chromosomes = df["chrom"].tolist()
    if int(subprocess.check_output(f"wc -l {bedFile} | cut -d ' ' -f 1", shell = True).decode("utf-8").rstrip()) != 0:
        chromosomes = set(pd.read_csv(bedFile, sep="\t", header=None)[0])
        chromosomes = [f"{chromosome}_region" for chromosome in chromosomes]
    df = df[df["chrom"].isin(chromosomes)]
    bases = df["bases"].sum()
    cov = df["mean"].mean()
    return [bases / 10**9, round(cov, 2)]

def parseMosDepthRegion(mosdepthFile, bedFile):
    if int(subprocess.check_output(f"wc -l {bedFile} | cut -d ' ' -f 1", shell = True).decode("utf-8").rstrip()) == 0:
        return ["NA"]
    if "offTarget" in bedFile:
        return [0]
    df = pd.read_csv(mosdepthFile, sep="\t", header=None, compression="gzip")
    try:
        df["depth"] = df[3].astype(float)
    except: #si bed est annote depth est dans la quatrieme colonne
        df["depth"] = df[4].astype(float)
    return [round(df[df["depth"] > 10].shape[0] / df.shape[0] * 100, 2)]

def computeQC(bamStatFile, summaryFile, regionFile, outFile, bed, wildcard):
    indexes = ["Longueur moyenne des reads (bp)", \
    "Nombre de bases sequencees (Gb)", \
    "Couverture moyenne", \
    "Pourcentage de région enrichie à plus de 10X"]
    dico = {}
    dico[wildcard] = parseBamStats(bamStatFile)
    dico[wildcard] += parseMosDepthSummary(summaryFile, bed)
    dico[wildcard] += parseMosDepthRegion(regionFile, bed)
    qc = pd.DataFrame(dico)
    qc.index = indexes
    qc.to_csv(outFile, sep="\t")

def concatQC(injections, run, outFile):
    injections = sorted(injections)
    final = pd.DataFrame()
    if len(injections) != 1:
        for injection in injections:
            final = pd.concat([final, pd.read_csv(injection, sep="\t", index_col=0)], axis=1)
    final = pd.concat([final, pd.read_csv(run, sep="\t", index_col=0)], axis=1)
    final.to_csv(outFile, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser(
	description="Compute QC for run")

    subs = parser.add_subparsers(dest = "command")
    computeQC_parser = subs.add_parser("qc", help="Compute qc from bamstat and coverage files")

    computeQC_parser.add_argument("-b", "--bed", help="Bed file.")
    computeQC_parser.add_argument("-B", "--bamStatFile", required=True, help="Bamstat file.")
    computeQC_parser.add_argument("-O", "--outFile", required=True, help="Path to output file.")
    computeQC_parser.add_argument("-r", "--mosdepthRegions", required=True, help="Mosdepth region file.")
    computeQC_parser.add_argument("-s", "--mosdepthSummary", required=True, help="Mosdepth summary file.")
    computeQC_parser.add_argument("-w", "--wildcard", required=True, help="Wildcard ID.")

    concat_parser = subs.add_parser("concat", help='Concate qc computed')
    concat_parser.add_argument("-i", "--injections", required=True, help="Path to qc injections files.", nargs="+")
    concat_parser.add_argument("-O", "--outFile", required=True, help="Path to output file.")
    concat_parser.add_argument("-r", "--run", required=True, help="Path to qc run file.")

    args = parser.parse_args()
    if args.command == "qc":
        computeQC(args.bamStatFile, args.mosdepthSummary, args.mosdepthRegions, args.outFile, args.bed, args.wildcard)
    if args.command == "concat":
        concatQC(args.injections, args.run, args.outFile)
