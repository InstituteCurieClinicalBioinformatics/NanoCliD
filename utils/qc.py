import subprocess
import pandas as pd
from argparse import ArgumentParser
import os
import re
import math

##pour calculer le N50, il suffit de faire un bedtools merge après avoir refait tourner le deuxième intersect.
##on créé une colonne avec la longueur des contigs.
##On somme cette longueur.
##On regarde pour quelle longueur de reads on arrive à somme / 2 -> N50. Pour ADAPTIVE_11 N50 = 226154

def parseBedtoolsIntersect(intersectFile):
    df = pd.read_csv(intersectFile, sep="\t", header=None)
    df[4] = df[2] - df[1]
    df = df.sort_values(by=[4], ascending=False)
    totalBases = df[4].sum()
    bases = 0
    for index, row in df.iterrows():
        bases += row[4]
        if bases >= totalBases / 2:
            n50 = row[4]
            break
    return [n50]
    # df[5] = df[3] * df[4]
    # return [math.trunc(df[5].sum() / df[3].sum())]

def parseMosDepthSummary(mosdepthFile, bedFile):
    chromosomes = set(pd.read_csv(bedFile, sep="\t", header=None)[0])
    chromosomes = [f"{chromosome}_region" for chromosome in chromosomes]
    df = pd.read_csv(mosdepthFile, sep="\t")
    df = df[df["chrom"].isin(chromosomes)]
    bases = df["bases"].sum()
    cov = df["mean"].mean()
    return [bases, cov]

def parseMosDepthRegion(mosdepthFile):
    df = pd.read_csv(mosdepthFile, sep="\t", header=None, compression="gzip")
    return [round(df[df[3] > 10].shape[0] / df.shape[0] * 100, 2)]

def computeQC(bedtoolsFile, summaryFile, regionFile, outFile, bed, wildcard):
    indexes = ["N50 (bp)", \
    "Nombre de bases sequencees", \
    "Couverture moyenne", 
    "Pourcentage de région enrichie à plus de 10X"]
    dico = {}
    dico[wildcard] = parseBedtoolsIntersect(bedtoolsFile)
    dico[wildcard] += parseMosDepthSummary(summaryFile, bed)
    dico[wildcard] += parseMosDepthRegion(regionFile)
    qc = pd.DataFrame(dico)
    qc.index = indexes
    qc.to_csv(outFile, sep="\t")

def concatQC(injections, run, outFile):
    injections = sorted(injections)
    final = pd.DataFrame()
    for injection in injections:
        final = pd.concat([final, pd.read_csv(injection, sep="\t", index_col=0)], axis=1)
    final = pd.concat([final, pd.read_csv(run, sep="\t", index_col=0)], axis=1)
    final.to_csv(outFile, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser(
	description="Compute QC for run")

    subs = parser.add_subparsers(dest = "command")
    computeQC_parser = subs.add_parser("qc", help="Compute qc from bamstat and coverage files")

    computeQC_parser.add_argument("-b", "--bed", required=True, help="Bed file.")
    computeQC_parser.add_argument("-i", "--intersectFile", required=True, help="Bedtools intersect file.")
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
        computeQC(args.intersectFile, args.mosdepthSummary, args.mosdepthRegions, args.outFile, args.bed, args.wildcard)
    if args.command == "concat":
        concatQC(args.injections, args.run, args.outFile)
