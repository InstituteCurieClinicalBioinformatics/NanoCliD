import pandas as pd
import subprocess
import os
from argparse import ArgumentParser

def openVCF(path):
    grep = "zgrep" if ".gz" in path else "grep"
    comments = subprocess.check_output(f"{grep} ^# {path}", shell=True).decode("utf-8").rstrip().split("\n")
    header = comments[-1].split("\t")
    df = pd.DataFrame(columns = header)
    compression = "gzip" if ".gz" in path else None
    if subprocess.call(f"{grep} -v '#' {path}", shell=True) != 1:
        df = pd.read_csv(path, sep = "\t", header = None, comment = "#", compression = compression)
        df.columns = header
    return comments[:-1], df

def cleanNanoCaller(path, outFile):
    comments, df = openVCF(path)
    outVcf = '.'.join(outFile.split('.')[:-1])
    if df.shape[0] != 0:
        df["ALT1"] = df["ALT"].str.split(",", expand=True)[0]
        if 1 in df["ALT"].str.split(",", expand=True).columns:
            df["ALT2"] = df["ALT"].str.split(",", expand=True)[1]
            df = df[~(df["ALT1"] == df["ALT2"])]
            df.pop("ALT2")
        df.pop("ALT1")
        f = open("tmp.vcf", "w")
        f.write("\n".join(comments))
        f.write("\n")
        f.close()
        df.to_csv("tmp2.vcf", index = None, sep = "\t")
        subprocess.call("cat tmp2.vcf >> tmp.vcf && rm -f tmp2.vcf", shell=True)
        subprocess.call(f"mv tmp.vcf {outVcf}", shell=True)
        subprocess.call(f"bgzip {outVcf} && tabix {outFile}", shell=True)
    else:
        f = open(outVcf, "w")
        f.write("\n".join(comments))
        f.write("\n")
        f.close()
        df.to_csv(outVcf, index = None, sep = "\t", mode = "a")
        subprocess.call(f"bgzip {outVcf} && tabix {outFile}", shell=True)
        subprocess.call(f"rm -f {outVcf}", shell = True)

if __name__ == "__main__":
    parser = ArgumentParser(description="Clean nanoCaller vcf. Remove lines when ALT is duplicated.")

    parser.add_argument("-i", "--inFile", required=True, help="Nanocaller vcf.")
    parser.add_argument("-o", "--outFile", required=True, help="Nanocaller cleaned vcf path.")

    args = parser.parse_args()

    cleanNanoCaller(args.inFile, args.outFile)
