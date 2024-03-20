import glob
import os
import subprocess
from argparse import ArgumentParser
import pandas as pd
import re
try:
   from utils.utils import getOutputs
except:
   pass
try:
    from utils.curie.utils import *
    curieFunctions = curieUtils()
    curieNetwork = True
except:
    curieNetwork = False

class NanoClid:

    def __init__(self, inputFolder=None, bedDir=None, run=None, outDir=None, dryRun=None, genomeVersion=None, until=None, samples="", outTemplate=None, snpEffDir=None, copyToTransverse=None, copyToWorkspace=None, sampleSheet = None, refDir = None, hostName = None, snakemakeBin = None, email = "", profile = None):
        self.inputFolder = inputFolder
        self.bedDir = bedDir
        self.run = run
        self.outDir = outDir
        self.dryRun = dryRun
        self.genomeVersion = genomeVersion
        self.until = until
        self.samples = samples.split(",")
        self.outTemplate = outTemplate
        self.gitDir = os.path.dirname(os.path.realpath(__file__))
        self.snpEffDir = snpEffDir
        if snakemakeBin is None and profile is None and curieNetwork:
            self.snakemakeBin, profile = curieFunctions._setSnakemakeBinAndProfile()
        else:
            self.snakemakeBin, profile = snakemakeBin, profile
        self.configTemplate = self.__setConfigTemplate(self.gitDir, self.run, curieNetwork, profile)
        self.profile = os.path.join(self.gitDir, "profiles", "curie", profile) if curieNetwork else os.path.join(self.gitDir, "profiles", "externe", profile)
        if refDir == "" and curieNetwork:
            self.refDir = curieFunctions._setRefDir(refDir, genomeVersion, profile)
        else:
            self.refDir = refDir
        self.copyToTransverse = copyToTransverse
        self.copyToWorkspace = copyToWorkspace
        self.sampleSheet = sampleSheet
        self.analysis = None
        self.combinaison = "''"
        self.fastqConcatenated = False
        self.fromBlow5 = "'no'"
        self.fromPod5 = "''"
        self.hostName = hostName if hostName else ""
        self.email = email
        if curieNetwork:
            self.email = "bioinfo-clinique@curie.fr"

    def __setConfigTemplate(self, gitDir, run, curieNetwork, profile):
        folder = "curie" if curieNetwork else "externe"
        if run == "A000":
            return os.path.join(gitDir, "data", folder, f"test_{profile}.yaml")
        return os.path.join(gitDir, "config", folder, profile, "config.yaml")

    def __checkBedIntegrity(self, bed):
        df = pd.read_csv(bed, sep = "\t", header = None)
        if not "chr" in df[0].iloc[0]:
            df = df.sort_values(by = [0,1])
            df[0] = df[0].astype(str)
            df[0] = "chr" + df[0]
            df.to_csv(bed, sep = "\t", index = None, header = False)

    def parseSampleSheet(self, sampleSheet):
        subprocess.call(f'mkdir -p {os.path.join(self.inputFolder, self.run, "archive")}', shell = True)
        if len(glob.glob(f"{self.inputFolder}/{self.run}/archive/*bed")) != 0:
            self.bed = glob.glob(f"{self.inputFolder}/{self.run}/archive/*bed")[0]
        else:
            code = subprocess.call(f'grep TargetBED {sampleSheet}', shell = True)
            if code == 1:
                subprocess.call(f"touch {self.inputFolder}/{self.run}/archive/empty.bed", shell = True)
                self.bed = f"{self.inputFolder}/{self.run}/archive/empty.bed"
            else:
                if os.path.basename(sampleSheet) == "A000_samplesheet.csv":
                    bed = os.path.join(self.gitDir, "data", "test.bed")
                else:
                    bed = subprocess.check_output(f'grep TargetBED {sampleSheet}', shell=True).decode('utf-8').rstrip().split(',')[1]
                    if curieNetwork:
                        bed = curieFunctions._getBed(bed)
                    else:
                        bed = os.path.join(self.bedDir, bed)
                subprocess.call(f"cp {bed} {self.inputFolder}/{self.run}/archive/", shell = True)
                self.bed = f"{self.inputFolder}/{self.run}/archive/{os.path.basename(bed)}"
                if not ".bed" in self.bed:
                    self.bed = f"{self.bed}.bed"
        if self.bed != f"{self.inputFolder}/{self.run}/archive/empty.bed":
            self.__checkBedIntegrity(self.bed)
        self.sequencer = ""
        if "M" in self.run:
            self.sequencer = "_mk1c"
        self.flowCellType = subprocess.check_output(f"grep FlowcellType {sampleSheet}", shell = True).decode("utf-8").rstrip().split(",")[1].lower()
        self.barcodingKits = subprocess.check_output(f"grep Assay {sampleSheet}", shell = True).decode("utf-8").rstrip().split(",")[1]
        skipRows = subprocess.check_output(f'grep -n "\[Data\]" {sampleSheet}', shell = True).decode("utf-8").split(":")[0]
        sampleSheet = pd.read_csv(sampleSheet, sep  = ",", skiprows = int(skipRows), keep_default_na = False)
        sampleSheet = sampleSheet[sampleSheet["Index_ID"] != ""]
        self.mergedSamples = sampleSheet["Index_ID"].tolist()
        sampleSheet["Index_ID"] = sampleSheet["Index_ID"].str[2:]
        sampleSheet["Index_ID"] = "barcode" + sampleSheet["Index_ID"]
        self.samplesPath = dict(zip(self.mergedSamples, sampleSheet.Index_ID))
        self.__setTypeAnalysis(self.inputFolder, self.run, self.mergedSamples)
        self.__setDataPath()

    def __setWildcardsCombinaison(self, injections, fast5Dir):
        self.combinaison = []
        injectionsForMerge = []
        samples = sorted(fast5Dir.keys())
        for injection in injections:
            injectionsForMerge.append(injection)
            for sample in samples:
                if injection in fast5Dir[sample] and "_" in sample:
                    self.combinaison.append({"samples" : sample, "injections" : injection})

    '''
    Sometimes, demultiplexing is not performed, so we need to perform demultiplexing
    '''
    def __isDemultiplexingToDo(self, fast5Paths):
        self.demultiplexing = False
        fast5 = []
        for fast5Pass in fast5Paths.values():
            if isinstance(fast5Pass, list):
                fast5.append(",".join(fast5Pass))
            else:
                fast5.append(fast5Pass)
        fast5 = set(",".join(fast5).split(",")) #get all fast5 paths
        if self.multiplexing and len(fast5) <= len(self.injections):
            self.demultiplexing = True

    '''
    Sometimes, basecalling can be performed on minknown. So if we have the fastq folder, just concat them in order to avoid to redo basecalling
    '''
    def __setFastqPath(self, folders):
        r = re.compile("fastq_pass")
        if len(list(filter(r.search, folders))) != 0:
            return list(filter(r.search, folders))[0]
        return ""

    def __concatFastq(self, fastqPaths, combinaison):
        for sample in fastqPaths.keys():
            if fastqPaths[sample] != "":
                resFolder = f"{self.outDir}"
                if self.analysis == "multipleInjections":
                    resFolder = os.path.join(f"{self.outDir}/{self.run}")
                injection = [dico["injections"] for dico in combinaison if dico["samples"] == sample][0]
                subprocess.call(f"mkdir -p {resFolder}", shell = True)
                subprocess.call(f"mkdir -p {resFolder}/{injection}", shell = True)
                subprocess.call(f"mkdir -p {resFolder}/{injection}/{sample}", shell = True)
                subprocess.call(f"mkdir -p {resFolder}/{injection}/{sample}/FASTQ", shell = True)
                subprocess.call(f"cat {fastqPaths[sample]}/*fastq.gz > {resFolder}/{injection}/{sample}/FASTQ/{sample}.fastq.gz", shell = True)
                self.fastqConcatenated = True

    def __setFast5Path(self, extension):
        self.fast5Paths = {}
        self.fastqPaths = {}
        self.samplesToMerge = {}
        self.samples = []
        samples = sorted(list(self.samplesPath.keys()))
        runDir = f"{self.inputFolder}" if self.analysis == "simpleInjection" else os.path.join(self.inputFolder ,self.run)
        if self.analysis == "simpleInjection":
            runDir = f"{self.inputFolder}"
        for sample in samples:
            if self.multiplexing or subprocess.call(f'find {runDir}/ -type d -regextype posix-egrep -regex ".*(pod5|fast5).*" | grep barcode', shell = True) == 0:
                pattern = self.samplesPath[sample]
            else:
                pattern = f"{extension}*"
            self.samplesToMerge[sample] = []
            self.fast5Paths[sample] = []
            for i in range(len(self.injections)):
                fast5Dir = self.getFiles(os.path.join(runDir, (self.injections[i])), 'd', pattern, True)
                fastqDir = [f for f in fast5Dir if "fastq" in f and not "fail" in f and not "skip" in f]
                fast5Dir = [f for f in fast5Dir if extension in f and not "fail" in f and not "skip" in f]
                if len(fast5Dir) == 0 and pattern == self.samplesPath[sample]: #if multiplexing and no fast5 with barcodeNb pattern found look for fast5 folder supposing demultiplexing is needed
                    fast5Dir = self.getFiles(os.path.join(runDir, (self.injections[i])), 'd', f'{extension}*', True)
                    fastqDir = [f for f in fast5Dir if "fastq" in f and not "fail" in f and not "skip" in f]
                    fast5Dir = [f for f in fast5Dir if extension in f and not "fail" in f and not "skip" in f]
                if len(fast5Dir) > 0:
                    self.fastqPaths[f"{sample}_{i+1}"] = self.__setFastqPath(fastqDir)
                    self.fast5Paths[sample].append(fast5Dir[0])
                    self.fast5Paths[f"{sample}_{i+1}"] = fast5Dir[0]
                    self.samples.append(f"{sample}_{i+1}")
                    self.samplesToMerge[sample].append(f"{sample}_{i+1}")
        samplesToRemove = []
        for sample in self.fast5Paths.keys():
            if len(self.fast5Paths[sample]) == 0:
                samplesToRemove.append(sample)
        for sample in samplesToRemove:
            del self.fast5Paths[sample]
            if not "_" in sample:
                del self.samplesToMerge[sample]
        self.mergedSamples = [sample for sample in self.fast5Paths.keys() if "_" not in sample]

    def __setReportFilesPath(self, fast5Paths, injections):
        self.reportFiles = {}
        for sample in fast5Paths.keys():
            if self.analysis == "simpleInjection":
                self.reportFiles[sample] = self.getFiles(os.path.join(self.inputFolder, self.run), "f", "report*.md", False)
            else:
                for injection in injections:
                    if injection in fast5Paths[sample]:
                        reportFile = self.getFiles(os.path.join(self.inputFolder, self.run, injection), "f", "report*.md", False)
                        if reportFile == "":
                            self.reportFiles[sample] = '""'
                        else:
                            self.reportFiles[sample] = reportFile

    def __setDataPath(self):
        self.__setFast5Path("pod5")
        self.fromPod5 = "'no'" if self.fast5Paths == {} else "'yes'"
        if self.fast5Paths == {} and self.fastqPaths == {}:
            self.__setFast5Path("fast5")
            if self.fast5Paths == {}:
                self.__setFast5Path("blow5")
                self.fromBlow5 = "'yes'"
        self.__setReportFilesPath(self.fast5Paths, self.injections)
        self.__setWildcardsCombinaison(sorted(self.injections), self.fast5Paths)
        self.__isDemultiplexingToDo(self.fast5Paths)
        self.__concatFastq(self.fastqPaths, self.combinaison)

    def __setTypeAnalysis(self, inputFolder, run, samples):
        self.multiplexing = False
        if len(samples) > 1:
            self.multiplexing = True
        folders = [folder for folder in os.listdir(f"{inputFolder}/{run}/") if os.path.isdir(f"{inputFolder}/{run}/{folder}")]
        injections = [folder for folder in folders if re.match(f"{run}_[0-9][0-9]*", folder)]
        self.analysis = "simpleInjection"
        if len(injections) > 0:
            self.analysis = "multipleInjections"
            self.injections = sorted([f"{injection}" for injection in injections])
        else:
            self.injections = [run]


    def parseLineConfig(self, line):
        if line.rstrip() != "":
            key, value = line.rstrip().split(":", 1)[0], line.rstrip().split(":", 1)[1]
            value = value.replace(" ", "", 1)
            if key[-1] == " ":
                key = key.replace(" ", "", 1)
            return key, value
        return "", ""

    def loadConfig(self, yaml):
        config = {}
        lines = open(yaml, "r").readlines()
        i = 0
        while i <= len(lines) - 1:
            key, value = self.parseLineConfig(lines[i])
            if value != "":
                config[key] = value
                i += 1
            elif key != "" and value == "":
                config[key] = {}
                i += 1
                while i < len(lines) and self.parseLineConfig(lines[i])[1] != "":
                    subKey, value = self.parseLineConfig(lines[i])
                    config[key][subKey.strip()] = value
                    i += 1
            else:
                i += 1
        return config

    def writeConfig(self, dico, path):
        f = open(path, "w")
        for key in dico.keys():
            if type(dico[key]) != dict:
                f.write(f"{key}: {dico[key]}\n")
            else:
                f.write(f'{key}: \n')
                for subKey in dico[key].keys():
                    f.write(f" {subKey}: {dico[key][subKey]}\n")
            f.write("\n")
        f.close()

    def __updateConfig(self, config):
        config["analysis"] = self.analysis
        config["demultiplexing"] = self.demultiplexing
        config["guppy"]["parameters"] = config["guppy"]["parameters"].replace("FLOWCELL", f"dna_{self.flowCellType}_e8.2_400bps_hac{self.sequencer}.cfg" if "r10" in self.flowCellType else f"dna_{self.flowCellType}_450bps_hac{self.sequencer}.cfg")
        if self.demultiplexing:
            config["guppy"]["parameters"] = f'{config["guppy"]["parameters"][1:-1]} --barcode_kits "{self.barcodingKits}"' #[1:-1] to remove quote
        config["flowCellType"] = self.flowCellType.split(".")[0]
        config["bed"] = self.bed
        config["sampleSheet"] = self.sampleSheet
        config["run"] = self.run
        config["genome_version"] = self.genomeVersion
        config["genome"] = os.path.join(self.refDir, f'{self.genomeVersion}.fasta') if os.path.exists(os.path.join(self.refDir, f'{self.genomeVersion}.fasta')) else os.path.join(self.refDir, f'{self.genomeVersion}.fa')
        config["genomeFile"] = os.path.join(self.refDir, f"{self.genomeVersion}.genome")
        config["minimap2"]["mmi"] = f"{self.refDir}/{self.genomeVersion}.mmi"
        config["wildcards"]["samples"] = self.samples
        config["wildcards"]["mergedSamples"] = self.mergedSamples
        config["wildcards"]["run"] = [self.run]
        config["wildcards"]["injections"] = self.injections
        config["samplesToMerge"] = self.samplesToMerge
        config["combinaison"] = self.combinaison
        config["git_dir"] = self.gitDir
        config["input_dir"] = self.inputFolder
        config["output_dir"] = f"{self.outDir}"
        config["template"] = self.outTemplate
        config["reportFiles"] = self.reportFiles
        config["fast5Dir"] = self.fast5Paths
        config["snpEff"]["dataDir"] = self.snpEffDir
        config["fromBlow5"] = self.fromBlow5
        config["fromPod5"] = self.fromPod5
        if self.email:
            config["email"] = self.email
        if curieNetwork:
            config = curieFunctions._addSpecificCurieInfoToConfig(config, self.run, self.hostName, self.gitDir)
        return config

    def __createConfig(self, template, folder, run):
        os.makedirs(folder, exist_ok=True)
        configFile = f"{run}.yaml"
        config = self.loadConfig(template)
        config = self.__updateConfig(config)
        self.writeConfig(config, os.path.join(folder, configFile))
        return os.path.join(folder, configFile)

    def __updateProfile(self, profile, outDir, run):
        profileType = profile.split("/")[-1]
        subprocess.call(f"cp -r {profile} {outDir}/{run}", shell=True)
        profileDico = self.loadConfig(f"{outDir}/{run}/{profileType}/config.yaml")
        profileDico["singularity-args"] = profileDico["singularity-args"].replace("GIT_DIR", self.gitDir)
        profileDico["singularity-args"] = profileDico["singularity-args"].replace("REF_DIR", self.refDir)
        profileDico["singularity-prefix"] = profileDico["singularity-prefix"].replace("GIT_DIR", self.gitDir)
        if "cluster" in profileDico.keys():
            profileDico["cluster"] = profileDico["cluster"].replace("logs_cluster", f"{self.outDir}/{self.run}/logs_cluster")
        self.writeConfig(profileDico, f"{outDir}/{run}/{profileType}/config.yaml")
        if "standalone" in profile and not self.until:
            #standalone and no until means you want to run the pipe in standalone configuration, so we remove run_on_cluster.txt from expected output
            newTemplate = os.path.join(self.outDir, self.run, os.path.basename(self.outTemplate))
            subprocess.call(f"mkdir -p {os.path.join(self.outDir, self.run)}", shell = True)
            subprocess.call(f"cp {self.outTemplate} {newTemplate}", shell = True)
            subprocess.call(f"sed -i '/run_on_cluster.txt/d' {newTemplate}", shell = True)
            self.outTemplate = newTemplate
        profile = f"{outDir}/{run}/{profileType}"
        return profile

    def runSnakemake(self, snakemakeBin, snakefile, configFile, profile, dryRun):
        cmd = f"{snakemakeBin} -s {snakefile} --profile {profile} --configfile {configFile} -d {os.path.join(self.outDir, self.run)}"
        demultiplexingRule = os.path.join(self.gitDir, "workflow/rules/demultiplexing.snk")
        if self.demultiplexing:
            cmdDemultiplexing = f"{snakemakeBin} -s {demultiplexingRule} --profile {profile} --configfile {configFile} -d {os.path.join(self.outDir, self.run)}"
            cmd = f"{cmdDemultiplexing} -t && {cmd}"
            if "standalone" in profile:
                cmd = cmdDemultiplexing
            if self.until == "guppy":
                cmd = f"{cmdDemultiplexing}"
        if self.until and "standalone" in profile:
            cmd = " ".join((cmd, f"-U {self.until}"))
        if dryRun:
            cmds = cmd.split(" && ")
            cmds = [f"{cmd} -n" for cmd in cmds]
            cmd = cmds[0]
            if len(cmds) > 1:
                cmd = " && ".join(cmds)
            print(cmd)
            subprocess.call(f"{cmd}", shell=True)
            return
        cmdNanoclidFinished = f"touch {self.outDir}/{self.run}/nanoclid_done.txt"
        cmd = " && ".join((cmd, cmdNanoclidFinished))
        print(cmd)
        subprocess.call(f"{cmd}", shell=True)

    def getFiles(self, path, kind, pattern, files=False):
        find = subprocess.check_output(f"find -L {path} -type {kind} -name '{pattern}'", shell=True).decode("utf-8").rstrip().split("\n")
        if len(find) == 0:
            print(f"No files found in path {path} with pattern {pattern}")
            return ""
        elif files:
            if pattern == "fast5*":
                return [file for file in find if not "fail" in file]
            return sorted(find)
        else:
            return find[0]

    def getSpecificPath(self, runDir, pattern, kind, injections):
        summaryFiles = self.getFiles(runDir, kind, pattern, files=True)
        summaryFilesPath = {}
        for injection in injections:
            for summaryFile in summaryFiles:
                if injection in summaryFile:
                    summaryFilesPath[injection] = summaryFile
                    break
        if pattern == "report*.md":
            injectionWithNoReport = set(injection).difference(set(summaryFilesPath.keys()))
            for injection in injectionWithNoReport:
                summaryFilesPath[injection] = ""
        if summaryFilesPath == {}:
            for injection in injection:
                summaryFilesPath[injection] = ""
        if pattern != "report*.md":
            self.injections = list(set(injections).intersection(set(summaryFilesPath.keys())))
        return summaryFilesPath

    def _runNanoClid(self):
        if not "standalone" in self.profile and curieNetwork:
            curieFunctions._deleteOldRun(self.inputFolder)
        if not self.sampleSheet and curieNetwork:
            self.sampleSheet = curieFunctions._getSampleSheet(self.run, self.inputFolder, self.loadConfig(self.configTemplate), self.loadConfig(os.path.join(self.profile, "config.yaml")), self.profile)
        self.parseSampleSheet(self.sampleSheet)
        configFile = self.__createConfig(self.configTemplate, f"{self.outDir}/{self.run}", self.run)
        profile = self.__updateProfile(self.profile, self.outDir, self.run)
        if curieNetwork:
            curieFunctions._recoverSymlink(profile, self.inputFolder, self.run)
            curieFunctions.runSnakemake(self.snakemakeBin, os.path.join(self.gitDir, "workflow/Snakefile"), configFile, profile, self.dryRun, self.outDir, self.run, self.gitDir, self.fastqConcatenated, self.until, self.demultiplexing, self.copyToTransverse, self.copyToWorkspace)
        else:
            self.runSnakemake(self.snakemakeBin, os.path.join(self.gitDir, "workflow/Snakefile"), configFile, profile, self.dryRun)

    def _install(self, singularityFolder):
        print(singularityFolder)
        subprocess.call(f"mkdir -p {singularityFolder}", shell = True)
        gitDir = os.path.dirname(os.path.realpath(__file__))
        imagesPath = "http://xfer.curie.fr/get/moe5ZZ8PGPD/images.tar.gz"
        nanovarPath = "http://xfer.curie.fr/get/Al1MrbPdC0C/nanovar.tar.gz"
        print("Downloading images ...")
        code = subprocess.call(f"wget -P {singularityFolder}/ {imagesPath}", shell = True)
        if code == 1:
            raise ValueError(f"Download images from {imagesPath} failed. Please retry or contact us.")
        print("Download images OK")
        print("Untar archive ...")
        code = subprocess.call(f"tar -xvf {singularityFolder}/images.tar.gz -C {singularityFolder} && mv {singularityFolder}/images/* {singularityFolder} && rm -rf {singularityFolder}/images {singularityFolder}/images.tar.gz", shell=True)
        if code == 1:
            raise ValueError(f"Untar images failed. Please check download integrity")
        print("Untar archive OK")
        print("Download nanovar archive ...")
        code = subprocess.call(f"wget -P {gitDir}/annotations/ {nanovarPath}", shell = True)
        if code == 1:
            raise ValueError(f"Download nanovar archive failed from {nanovarPath}. Please retry and contact us")
        print("Download nanovar archive OK")
        print("Untar nanovar folder ...")
        code = subprocess.call(f"tar -xvf {gitDir}/annotations/nanovar.tar.gz -C {gitDir}/annotations/ && rm {gitDir}/annotations/nanovar.tar.gz", shell = True)
        if code == 1:
            raise ValueError(f"Untar nanovar folder failed. Please check git clone integrity. You should get git-lfs")
        print("Untar nanovar folder OK")
        print("Setting python virtual env ...")
        code1 = subprocess.call(f"pip3 install virtualenv", shell=True)
        code2 = subprocess.call(f"mkdir -p {gitDir}/venv && python3 -m venv {gitDir}/venv", shell=True)
        code3 = subprocess.run(['bash', '-c', f"source {gitDir}/venv/bin/activate && pip3 install --upgrade pip wheel setuptools"], check = True)
        code4 = subprocess.run(['bash', '-c', f"source {gitDir}/venv/bin/activate && pip3 install snakemake pandas"], check = True)
        if 1 in [code1, code2, code3, code4]:
            raise ValueError("Setting python virtualenv FAILED")
        print("Setting python virtualenv OK")
        print("Retrieving snpEff annotations ...")
        code = subprocess.call(f"wget -P {gitDir}/annotations/ http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19.zip && cd {gitDir}/annotations/ && unzip snpEff_v4_3_hg19.zip && rm snpEff_v4_3_hg19.zip", shell=True)
        if code == 1:
            raise ValueError("snpEff annotations retrieving failed for path http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19.zip. Link may have been removed. Please contact us.")
        print("Retrieving snpEff annotations OK")
        print("Install completed. NanoCliD is ready to use. Please launch test command.")

    def _runTest(self, refDir, profile, email):
        gitDir = os.path.dirname(os.path.realpath(__file__))
        print("Extracting data test folder...")
        code = subprocess.call(f"tar -xvf {os.path.join(gitDir, 'data', 'input', 'input.tar.gz')} -C {os.path.join(gitDir, 'data')}", shell = True)
        if code == 1:
            raise ValueError("Extracting data test folder OK")
        inputFolder = os.path.join(gitDir, "data")
        bedDir = os.path.join(gitDir, "data")
        run = "A000"
        outDir = os.path.join(gitDir, "data")
        dryRun = False
        genomeVersion = "hg19"
        until = False
        samples = ""
        outTemplate = os.path.join(gitDir, "data", "externe", "TEST.template")
        snpEffDir = os.path.join(gitDir, "annotations", "data")
        copyToTransverse = False
        copyToWorkspace = False
        sampleSheet = os.path.join(gitDir, "data", "A000_samplesheet.csv")
        hostName = None
        snakemakeBin = os.path.join(gitDir, "venv", "bin", "snakemake")
        nanoclid = NanoClid(inputFolder, bedDir, run, outDir, dryRun, genomeVersion, until, samples, outTemplate, snpEffDir, copyToTransverse, copyToWorkspace, sampleSheet, refDir, hostName, snakemakeBin, email, profile)
        nanoclid._runNanoClid()

if __name__ == "__main__":
    parser = ArgumentParser(
    description="Launch NanoClid pipeline")

    subs = parser.add_subparsers(dest = "command")

    install_parser = subs.add_parser("install", help='Install NanoClid images.')
    install_parser.add_argument("-S", "--singularityFolder", help="Where to download singularity images ? Default is GIT_DIR/singularity", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'singularity'))

    test_parser = subs.add_parser("test", help="Run test")
    test_parser.add_argument("-e", "--email", help="Email adress to send run informations.")
    test_parser.add_argument("-p", "--profile", required=True, help="Profile to use to launch NanoCliD. Must be standalone|calcsub|abacus")
    test_parser.add_argument("-R", "--refDir", required=True, help="Path to folder containing genome files")

    run_parser = subs.add_parser("run", help='Run NanoClid')
    run_parser.add_argument("-b", "--bedDir", help="Path to folder containing bed files.", default = "")
    run_parser.add_argument("-B", "--snakemakeBin", help="Path to snakemake bin.")
    run_parser.add_argument("-c", "--noCopyToTransverse", action='store_false', help="Do not copy results to transverse. Only for curie network.")
    run_parser.add_argument("-D", "--snpEffDir", help="Folder containing snpEff annotation files.",  default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'annotations/data'))
    run_parser.add_argument("-e", "--email", help="Used by snakemake to send error/launching message.")
    run_parser.add_argument("-g", "--genomeVersion", required=True, help="Genome version for the analysis.")
    run_parser.add_argument("--hostName", help="Cluster on which launch the analysis. Only for curie network")
    run_parser.add_argument("-I", "--inputFolder", required=True, help="Path to the folder containing nanopore sequencing data.")
    run_parser.add_argument("-n", "--dryRun", action='store_true', help="Run nanoclid in dry run mode.")
    run_parser.add_argument("-O", "--outDir", help="Path to output folder.")
    run_parser.add_argument("-p", "--profile", help="Profile to use to launch NanoCliD. Must be standalone|calcsub|abacus")
    run_parser.add_argument("-r", "--runID", required = True, help="runID")
    run_parser.add_argument("-R", "--refDir", help="Ref dir containing genome annotations and reference files.", default="")
    run_parser.add_argument("-s", "--samples", help="Only run analysis on these samples. Must be comma separated.", default="")
    run_parser.add_argument("-S", "--samplesheet", help="Path to the samplesheet.")
    run_parser.add_argument("-T", "--outputTemplate", help="Path to the output template.", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', 'externe', 'ADAPTIVE.template'))
    run_parser.add_argument("-U", "--until", help="Specify until which rule you want to run the workflow", default="")
    run_parser.add_argument("-w", "--noCopyToWorkspace", action='store_false', help="Do not copy results to workspace. Only for curie network")

    args = parser.parse_args()

    if args.command == "install":
        NanoClid(profile = "standalone")._install(args.singularityFolder)

    if args.command == "test":
        NanoClid(profile = "standalone")._runTest(args.refDir, args.profile, args.email)

    if args.command == "run":
        if not args.outDir:
            args.outDir = args.inputFolder
        nanoclid = NanoClid(args.inputFolder, args.bedDir, args.runID, args.outDir, args.dryRun, args.genomeVersion, args.until, args.samples, args.outputTemplate, args.snpEffDir, args.noCopyToTransverse, args.noCopyToWorkspace, args.samplesheet, args.refDir, args.hostName, args.snakemakeBin, args.email, args.profile)
        nanoclid._runNanoClid()
