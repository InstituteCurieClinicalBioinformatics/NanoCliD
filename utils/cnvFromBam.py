from argparse import ArgumentParser
from pathlib import Path
from cnv_from_bam import iterate_bam_file

def main(bamFile, threads, mapQ, widthFile):
    result = iterate_bam_file(Path(bamFile), threads=threads, mapq_filter=threads)
    with open(widthFile, "w") as binFile:
        binFile.write(str(result.bin_width))

if __name__ == "__main__":
    parser = ArgumentParser(description="CNV from BAM")

    parser.add_argument("-b", "--bam", required = True, help = "Bam file")
    parser.add_argument("-m", "--mapQ", type = int, default = 60)
    parser.add_argument("-t", "--threads", type = int, default = 4)
    parser.add_argument("-w", "--widthFile", required = True)

    args = parser.parse_args()

    main(args.bam, args.threads, args.mapQ, args.widthFile)