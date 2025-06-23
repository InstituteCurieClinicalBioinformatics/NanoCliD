import pandas as pd
import os
from argparse import ArgumentParser

def extractTag(tagName, text, sortDescending=False):
    startIdxs = []
    endIdxs = []
    for i in range(len(text)):
        if f'<{tagName}' in text[i]:
            startIdxs.append(i)
        if f'</{tagName}' in text[i]:
            endIdxs.append(i)
    if sortDescending:
        startIdxs.sort(reverse=True)
        endIdxs.sort(reverse=True)
    return dict(zip(startIdxs, endIdxs))

def parseTable(idxs, text, qcLines):
    for key in idxs.keys():
        subText = text[key:idxs[key] + 1]
        idxsTr = extractTag("tr", subText, sortDescending=True) #pour recuperer la derniere entree tr dans la table plus facile de sort en descending
        startIdx, endIdx = list(idxsTr.keys())[0], list(idxsTr.values())[0]
        subText.insert(endIdx+1, qcLines)
        text[key:idxs[key] + 1] = subText
    return text

def convertQcToHtml(qcFile, sample):
    qcType = os.path.basename(qcFile).split("_qc_")[-1].split(".")[0]
    df = pd.read_csv(qcFile, sep = "\t", index_col = 0)
    if qcType == "ontarget":
        df = df.drop([df.index[0]])
    linesTableHtml = ""
    for index, row in df.iterrows():
        linesTableHtml += f'''<tr>
    <th>{index} {qcType}</th>
    <td>{row[sample]}</td>
</tr>\n'''
    return linesTableHtml
        

def parseReport(htmlFile, tagName, qcFiles, sample):
    with open(htmlFile, "r") as f:
        lines = f.readlines()
    qcLinesHtml = ""
    for i in range(len(qcFiles)):
        qcLinesHtml += convertQcToHtml(qcFiles[i], sample)
    idxs = extractTag(tagName, lines)
    text = parseTable(idxs, lines, qcLinesHtml)
    with open(f"{htmlFile.replace('.tmp', '')}", "w") as f:
        f.write("".join(text))

if __name__ == "__main__":
    parser = ArgumentParser(description="Add qc table to NanoPlot report")

    parser.add_argument("-i", "--inFile", required=True, help="NanoPlot report.")
    parser.add_argument("-q", "--qcFiles", required=True, help="QC file from NanoCliD.", nargs = "+")
    parser.add_argument("-s", "--sampleID", required=True, help="Sample ID.")
    parser.add_argument("-t", "--tagName", help="Tag name of the table in the html", default = "table")

    args = parser.parse_args()

    parseReport(args.inFile, args.tagName, args.qcFiles, args.sampleID)