import pandas as pd
import subprocess
import os
from argparse import ArgumentParser

def concatSV(vcfs, outFile):
	final = pd.DataFrame()
	comments = []
	callers = []
	for vcf in vcfs:
		if subprocess.call(f"grep -cv '#' {vcf}", shell = True) != 1:
			header = subprocess.check_output(f"grep '^##' {vcf}", shell = True).decode("utf-8").rstrip().split("\t")
			comments += list(set(header).difference(set(comments)))
			columns = subprocess.check_output(f"grep '#CHROM' {vcf}", shell = True).decode("utf-8").rstrip().split("\t")
			df = pd.read_csv(vcf, sep = "\t", header = None, comment = "#")
			df.columns = columns
			final = pd.concat([df, final])
			callers += [os.path.dirname(vcf).split("/")[-1]] * df.shape[0]
	comments.append("\n")
	with open(outFile, "w") as fh:
		fh.write("".join(comments))
	final["Caller"] = callers
	final.to_csv(outFile, sep = "\t", index = None, mode = "a")


if __name__ == "__main__":
	parser = ArgumentParser(description="Concat SV vcfs")

	parser.add_argument("-v", "--vcfs", required = True, help = "VCFs files", nargs = "+")
	parser.add_argument("-o", "--outFile", required = True)

	args = parser.parse_args()

	concatSV(args.vcfs, args.outFile)
