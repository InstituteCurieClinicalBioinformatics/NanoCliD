import pandas as pd
from collections import defaultdict
import pycircos
import matplotlib.pyplot as plt
import collections
import os
from argparse import ArgumentParser

def parseAnnotated(path, gene, run, outDir):
    fusionColumns = ["chr1", "start1", "end1", "chr2", "start2", "end2"]
    columns = [f'sniffles.{run}', f'cuteSV.{run}', f'SVIM.{run}', f'nanovar.{run}']
    df = pd.read_csv(path, sep="\t")
    columns = list(set(set(columns).intersection(df.columns)))
    df = df[(df["SV_type"] == "TRA") & (df["Gene_name"] == gene)][columns]
    fusion = pd.DataFrame(columns = fusionColumns)
    for index, row in df.iterrows():
        variants = row[columns].tolist()
        coords = [variant.split(":")[-1] for variant in variants if variant.split(":")[-1] != "NAN"]
        dico = defaultdict(list)
        for coord in coords:
            if "," in coord:
                coords = coord.split(",")
                for coord in coords:
                    dico = parseCoord(coord, dico)
            else:
                dico = parseCoord(coord, dico)
            fusion = pd.concat([pd.DataFrame(data = [writeCoords(dico)], columns = fusionColumns), fusion])
    fusion.to_csv(f"{outDir}/{run}_{gene}.csv", sep=",", index=None)

def parseCoord(coord, dico):
    start, end = coord.split("-")
    dico[start.split("_")[0]].append(start.split("_")[1])
    dico[end.split("_")[0]].append(end.split("_")[1])
    return dico

def writeCoords(dico):
    chromA = list(dico.keys())[0]
    chromB = list(dico.keys())[1]
    if len(set(dico[chromA])) == 1:
        startA = int(dico[chromA][0])
        endA = int(dico[chromA][0]) + 1
        startB = int(dico[chromB][0])
        endB = int(dico[chromB][0]) + 1
    else:
        startA = int(sorted(dico[chromA])[0])
        endA = int(sorted(dico[chromA])[-1])
        startB = int(sorted(dico[chromB])[0])
        endB = int(sorted(dico[chromB])[-1])
    return [chromA, startA, endA, chromB, startB, endB]

def crossGTF(chr, start, end, gtf):
    gtf = gtf[(gtf[0] == chr) & (start >= gtf[1]) & (end <= gtf[2])]
    gene = list(set(gtf[3].str.extract(r'gene_name "(\w+)"')[0]))
    if len(gene) == 0:
        gene = ["NAN"]
    return gene

def setCircosPlot(chromosomeFile, cytobandeFile):
    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle

    circle = Gcircle()

    with open(chromosomeFile) as f:
        f.readline()
        for line in f:
            line   = line.rstrip().split(",") 
            name   = line[0]
            length = int(line[-1])
            arc    = Garc(arc_id=name, size=length, interspace=3, raxis_range=(850,900), labelposition=50, label_visible=True)
            circle.add_garc(arc)

    circle.set_garcs()

    color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
                "acen":"#D82322"}

    arcdata_dict = collections.defaultdict(dict)

    with open(cytobandeFile) as f:
        f.readline()
        for line in f:
            line  = line.rstrip().split(",")
            name  = line[0]
            start = int(line[1])-1 
            width = int(line[2])-(int(line[1])-1) 
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["widths"]    = [] 
                arcdata_dict[name]["colors"]    = [] 
            arcdata_dict[name]["positions"].append(start) 
            arcdata_dict[name]["widths"].append(width)
            arcdata_dict[name]["colors"].append(color_dict[line[-1]])


    for key in arcdata_dict:
        circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"], width=arcdata_dict[key]["widths"], raxis_range=[850,900], facecolor=arcdata_dict[key]["colors"]) 

    return circle

def makeCircosPlot(chromosomeFile, cytobandeFile, fusionFile, run, outDir, gene, multipleFusions, gtf):

    circle = setCircosPlot(chromosomeFile, cytobandeFile)

    fusions = pd.read_csv(fusionFile, sep=",")

    legend = []

    for index, row in fusions.iterrows():
        chr1 = row["chr1"]
        start1 = row["start1"]
        end1 = row["end1"]
        chr2 = row["chr2"]
        start2 = row["start2"]
        end2 = row["end2"]
        gene1 = crossGTF(chr1, start1, end1, gtf)
        gene2 = crossGTF(chr2, start2, end2, gtf)
        source = (chr1, start1, end1, 900)
        destination = (chr2, start2, end2, 900)
        facecolor=circle.garc_dict[chr2].facecolor
        circle.chord_plot(source, destination, facecolor=facecolor, edgecolor=facecolor, linewidth=1)
        legend.append(f"{'|'.join(gene1)}-{'|'.join(gene2)}")
        if not multipleFusions:
            circle.figure.set_size_inches(8,9)
            circle.figure.suptitle(f"{run} {'|'.join(gene1)}-{'|'.join(gene2)}", y=0.97)
            plt.legend(legend, loc='lower center', bbox_to_anchor=(0.5, -0.05))
            if fusions.shape[0] == 1:
                save = f"{outDir}/{run}_{'|'.join(gene1)}-{'|'.join(gene2)}.png"
            else:
                save = f"{outDir}/{run}_{'|'.join(gene1)}-{'|'.join(gene2)}_{index}.png"
            circle.figure.savefig(save, dpi=300)
            circle = setCircosPlot(chromosomeFile, cytobandeFile)
            legend = []
    if multipleFusions:
        circle.figure.set_size_inches(8,10)
        circle.figure.suptitle(f"{run} {'|'.join(gene1)}-{'|'.join(gene2)}", y=0.97)
        plt.legend(legend, loc='lower center', bbox_to_anchor=(0.5, -0.05))
        circle.figure.savefig(f"{outDir}/{run}_{'|'.join(gene1)}-{'|'.join(gene2)}.png")
        circle = setCircosPlot(chromosomeFile, cytobandeFile)

def main(annotatedFile, genes, run, outDir, chromosomeFile, cytobandeFile, multipleFusions, gtf):
    if os.path.isfile(genes):
        genes = pd.read_csv(genes, header=None)[0].tolist()
    else:
        genes = gene.split(",")
    gtf = pd.read_csv(gtf, sep="\t", comment="#", header=None)
    gtf = gtf[[0,2,3,4,8]]
    gtf = gtf[gtf[2] == "gene"]
    gtf.pop(2)
    gtf.columns = [i for i in range(len(gtf.columns))]
    for gene in genes:
        parseAnnotated(annotatedFile, gene, run, outDir)
        makeCircosPlot(chromosomeFile, cytobandeFile, f"{outDir}/{run}_{gene}.csv", run, outDir, gene, multipleFusions, gtf)

if __name__ == "__main__":
    parser = ArgumentParser(description="Make circos plot")

    parser.add_argument("-b", "--cytobandeFile", required=True, help="File containing cytobande coordinates.")
    parser.add_argument("-c", "--chromosomeFile", required=True, help="File containing chromosomes length.")
    parser.add_argument("-g", "--gene", help="Gene to looking for genomic event.")
    parser.add_argument("-G", "--gtf", required=True, help="GTF file to know which genes are involved in genomic event.")
    parser.add_argument("-i", "--annotatedFile", required=True, help="Path to annotated SV file.")
    parser.add_argument("-l", "--listGenes", help="File containing genes to look for genomic event.")
    parser.add_argument("-m", "--multipleFusions", help="If one gene is related to multiple genomic events, they will all be printed in the same circos. Default is one circos per genomic event.", action='store_true')
    parser.add_argument("-O", "--outputDir", required=True, help="Path to output folder.")
    parser.add_argument("-r", "--run", required=True, help="Run ID.")

    args = parser.parse_args()


    if not args.gene and not args.listGenes:
        print("ERROR, you must use gene arguments or listGenes arguments.")
    else:
        if args.gene:
            main(args.annotatedFile, args.gene, args.run, args.outputDir, args.chromosomeFile, args.cytobandeFile, args.multipleFusions, args.gtf)
        else:
            main(args.annotatedFile, args.listGenes, args.run, args.outputDir, args.chromosomeFile, args.cytobandeFile, args.multipleFusions, args.gtf)