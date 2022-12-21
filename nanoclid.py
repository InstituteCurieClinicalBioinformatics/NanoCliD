import glob
import os
import subprocess
from argparse import ArgumentParser

class NanoClid:
    _SNAKEMAKE_MINION_PARAMETERS = "-p --verbose --latency-wait 60 --keep-going --use-singularity"
    _SINGULARITY_TAR_PART1 = "1rnOcr8M_lH3KSu3dIKpi6C06cxSt3zMN"
    _SINGULARITY_TAR_PART2 = "1-ZkxulbbVv54y__Zz6y-Xahh-YwftPUm"
    _ANNOTATIONS_TAR = "1RQa5QHgMDmcFR-IRSKwbUiMn_QMMoYKc"
    _INPUT_TAR = "1LpguQ1aUFQPGcnh4iow73rExaHT-E96_"
    _EXPECTED_TAR = "1hkn4YxjmID30Rqil58XjdiEUzm_1RTLq"
    _CMD_DOWNLOAD = """wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=FILEID" -O FILENAME && rm -rf /tmp/cookies.txt"""

    def __init__(self, inputFolder=None, run=None, outDir=None, bed=None, dryRun=None, genomeVersion=None, gtf=None, geneList=None, snakemakeParameters=None, cores=None, bind=None, refDir=None, samples="", configTemplate=None, outTemplate=None, gitDir=None, containersPath=None, snpEffDir=None, cpu=None):
        self.inputFolder = inputFolder
        self.run = run
        self.outDir = outDir
        self.bed = bed
        self.dryRun = dryRun
        self.genomeVersion = genomeVersion
        self.gtf = gtf
        self.geneList = geneList
        self.snakemakeParameters = snakemakeParameters
        self.cores = cores
        self.bind = bind
        self.refDir = refDir
        self.samples = samples.split(",")
        self.configTemplate = configTemplate
        self.outTemplate = outTemplate
        self.gitDir = gitDir
        self.containersPath = containersPath
        self.snpEffDir = snpEffDir
        self.cpu = cpu

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
        while i < len(lines) - 1:
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
                f.write("\n")
                f.write(f'{key}: \n')
                for subKey in dico[key].keys():
                    f.write(f" {subKey}: {dico[key][subKey]}\n")
        f.close()

    def __updateConfig(self, config):
        if self.cpu:
            config["guppy"]["runner"] = config["guppy"]["runner_cpu"]
            config["guppy"]["parameters"] = config["guppy"]["parameters_cpu"]
        else:
            config["guppy"]["runner"] = config["guppy"]["runner_gpu"]
            config["guppy"]["parameters"] = config["guppy"]["parameters_gpu"]
        config["bed"] = self.bed
        config["run"] = self.run
        config["genome_version"] = self.genomeVersion
        config["genome"] = f"{self.refDir}/{self.genomeVersion}.fasta"
        config["genomeFile"] = f"{self.refDir}/{self.genomeVersion}.genome"
        config["minimap2"]["mmi"] = f"{self.refDir}/{self.genomeVersion}.mmi"
        config["samples"] = self.samples
        config["git_dir"] = self.gitDir
        config["input_dir"] = self.inputFolder
        config["output_dir"] = self.outDir
        config["analysis"] = self.outTemplate
        config["containers_path"] = f"{self.containersPath}/"
        config["summaryFiles"] = self.summaryFiles
        config["reportFiles"] = self.reportFiles
        config["fast5Dir"] = self.fast5Dir
        config["circos"]["cytobande"] = f"{self.gitDir}/annotations/cytobande.csv"
        config["circos"]["chromosome"] = f"{self.gitDir}/annotations/chromosomes.tsv"
        config["circos"]["gtf"] = self.gtf
        config["circos"]["geneList"] = self.geneList
        config["snpEff"]["dataDir"] = self.snpEffDir
        return config

    def __createConfig(self, template, folder, run):
        os.makedirs(folder, exist_ok=True)
        configFile = f"{run}.yaml"
        config = self.loadConfig(template)
        config = self.__updateConfig(config)
        self.writeConfig(config, os.path.join(folder, configFile))
        return os.path.join(folder, configFile)

    def __checkDownloads(self, gitDir, files=False):
        defFiles = glob.glob(os.path.join(gitDir, "singularity", "*.def"))
        sifFiles = glob.glob(os.path.join(gitDir, "singularity", "*.sif"))
        defFiles = [os.path.basename(defFile).split(".")[0] for defFile in defFiles]
        sifFiles = [os.path.basename(sifFile).split(".")[0] for sifFile in sifFiles]
        if sorted(defFiles) != sorted(sifFiles):
            raise ValueError(f"Image for {set(defFiles).difference(sifFiles)} is missing.\nDownload must do not have ended properly.")

        annotations = ["hg19.dict", "hg19.fasta", "hg19.fasta.fai", "hg19.mmi", "hg19.genome", "hg19_subsampling.fa"]
        for annotation in annotations:
            if not os.path.isfile(os.path.join(gitDir, "annotations", annotation)):
                raise ValueError(f"{os.path.join(gitDir, 'annotations', annotation)} is missing\nDownload must do not have ended properly.")

        if files:
            files = open(os.path.join(gitDir, "data/input.template"), "r").readlines()
            files = [f.replace("INPUT", os.path.join(gitDir, "data/input")) for f in files]
            for f in files:
                if not os.path.isfile(f.rstrip()):
                    raise ValueError(f"{f} is missing\nDownload for input files must do not have ended properly.")

            files = open(os.path.join(gitDir, "data/check.template"), "r").readlines()
            filesCPU = [f.replace("FOLDER", os.path.join(gitDir, "data/expected/cpu")) for f in files]
            filesGPU = [f.replace("FOLDER", os.path.join(gitDir, "data/expected/gpu")) for f in files]
            files = filesCPU + filesGPU
            for f in files:
                if not os.path.isfile(f.rstrip()):
                    raise ValueError(f"{f} is missing\nDownload for expected files must do not have ended properly.")

    def runSnakemake(self, snakefile, configFile, arguments, bind, cores, dryRun):
        cmd = f"snakemake -s {snakefile} --configfile {configFile} {arguments} --singularity-args '-B {bind},{self.gitDir}/annotations/nanovar:/bnd --no-home --cleanenv --writable-tmpfs' -c {cores}"
        if dryRun:
            cmd = " ".join((cmd, "-n"))
        subprocess.call(f"{cmd}", shell=True)

    def getFiles(self, path, kind, pattern, files=False):
        find = subprocess.check_output(f"find {path} -type {kind} -name {pattern}", shell=True).decode("utf-8").rstrip().split("\n")
        if len(find) == 0:
            print(f"No files found in path {path} with pattern {pattern}")
            return ""
        elif files:
            if pattern == "fast5*":
                return [file for file in find if not "fail" in file]
            return find
        else:
            return find[0]

    def getSpecificPath(self, runDir, pattern, kind, samples):
        summaryFiles = self.getFiles(runDir, kind, pattern, files=True)
        summaryFilesPath = {}
        for sample in samples:
            for summaryFile in summaryFiles:
                if sample in summaryFile:
                    summaryFilesPath[sample] = summaryFile
                    break
        self.samples = list(set(samples).intersection(set(summaryFilesPath.keys())))
        return summaryFilesPath


    def _install(self, folder):
        gitDir = os.path.dirname(os.path.realpath(__file__))
        singularityDir = os.path.join(gitDir, 'singularity')
        defFiles = glob.glob(f"{singularityDir}/*.def")
        defFiles = sorted(defFiles)
        downloadDir = "/tmp/"
        if not folder:
            print("Downloading archives ...")
            subprocess.call(NanoClid._CMD_DOWNLOAD.replace("FILEID", NanoClid._SINGULARITY_TAR_PART1).replace("FILENAME", f"{downloadDir}/singularity1.tar.gz"), shell = True)
            subprocess.call(NanoClid._CMD_DOWNLOAD.replace("FILEID", NanoClid._SINGULARITY_TAR_PART2).replace("FILENAME", f"{downloadDir}/singularity2.tar.gz"), shell = True)
            subprocess.call(NanoClid._CMD_DOWNLOAD.replace("FILEID", NanoClid._ANNOTATIONS_TAR).replace("FILENAME", f"{downloadDir}/annotations.tar.gz"), shell = True)
            subprocess.call(NanoClid._CMD_DOWNLOAD.replace("FILEID", NanoClid._INPUT_TAR).replace("FILENAME", f"{downloadDir}/input.tar.gz"), shell = True)
            subprocess.call(NanoClid._CMD_DOWNLOAD.replace("FILEID", NanoClid._EXPECTED_TAR).replace("FILENAME", f"{downloadDir}/expected.tar.gz"), shell = True)
        else:
            downloadDir = folder
        print("Installation is running ...")
        subprocess.call(f"tar -xvf {downloadDir}/singularity1.tar.gz -C {downloadDir}", shell=True)
        subprocess.call(f"tar -xvf {downloadDir}/singularity2.tar.gz -C {downloadDir}", shell=True)
        subprocess.call(f"tar -xvf {downloadDir}/annotations.tar.gz -C {downloadDir}", shell=True)
        subprocess.call(f"tar -xvf {downloadDir}/input.tar.gz -C {os.path.join(gitDir, 'data')}/", shell=True)
        subprocess.call(f"tar -xvf {downloadDir}/expected.tar.gz -C {os.path.join(gitDir, 'data')}/", shell=True)
        subprocess.call(f"cp {downloadDir}/singularity*/*sif {os.path.join(gitDir, 'singularity')}/", shell = True)
        subprocess.call(f"cp {downloadDir}/annotations/* {os.path.join(gitDir, 'annotations')}/", shell = True)
        self.__checkDownloads(gitDir)
        subprocess.call(f"rm -rf {downloadDir}/singularity*.tar.gz && rm -rf {downloadDir}/singularity*/", shell = True)
        subprocess.call(f"rm -rf {downloadDir}/annotations.tar.gz && rm -rf {downloadDir}/annotations/", shell = True)
        subprocess.call(f"rm -rf {downloadDir}/input.tar.gz", shell = True)
        subprocess.call(f"rm -rf {downloadDir}/expected.tar.gz", shell = True)
        print("Downloading archives OK")
        print("Setting python virtualenv ...")
        code1 = subprocess.call(f"pip3 install virtualenv", shell=True)
        code2 = subprocess.call(f"mkdir -p {gitDir}/venv && python3 -m venv {gitDir}/venv", shell=True)
        code3 = subprocess.run(['bash', '-c', f"source {gitDir}/venv/bin/activate && pip3 install snakemake pandas python-circos"])
        if 1 in [code1, code2, code3]:
            print("Setting python virtualenv FAILED")
            return
        else:
            print("Setting python virtualenv OK")
            print("Unpackaging test data and annotations ...")
            codeTest = subprocess.call(f"tar -xvf {gitDir}/data/input/input.tar.gz -C {gitDir}/data/input/", shell=True)
            subprocess.call(f"mkdir -p {gitDir}/data/expected/cpu {gitDir}/data/expected/gpu", shell = True)
            codeExpected_cpu = subprocess.call(f"tar -xvf {gitDir}/data/expected/expected_cpu.tar.gz -C {gitDir}/data/expected/cpu", shell=True)
            codeExpected_gpu = subprocess.call(f"tar -xvf {gitDir}/data/expected/expected_gpu.tar.gz -C {gitDir}/data/expected/gpu", shell=True)
            self.__checkDownloads(gitDir, True)
            codeAnnotations = subprocess.call(f"tar -xvf {gitDir}/annotations/nanovar.tar.gz -C {gitDir}/annotations", shell=True)
            if 1 in [codeTest, codeExpected_cpu, codeExpected_gpu, codeAnnotations]:
                print("Unpackaging data FAILED")
                return
            else:
                print("Unpackaging data OK")
                print("Retrieving snpEff annotations ...")
                code = subprocess.call(f"cd {gitDir}/annotations && wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19.zip && unzip snpEff_v4_3_hg19.zip", shell=True)
                if code == 1:
                    print("Retrieving snpEff annotations FAILED")
                    return
                else:
                    if not os.path.isfile(os.path.join(f"{gitDir}/annotations/gencode.v40lift37.annotation.gtf.gz")):
                        print("Retrieving snpEff annotations OK")
                        print("Retrieving gtf file ...")
                        code = subprocess.call(f"cd {gitDir}/annotations && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz && gunzip gencode.v40lift37.annotation.gtf.gz", shell=True)
                        if code == 1:
                            print("Retrieving gtf file FAILED")
                            return
                        else:
                            print("Retrieving gtf file OK")

    def _runTest(self, cpu):
        print("Test is running ...")
        print("Checking if guppy_basecall_client is in PATH ...")
        if subprocess.call("which guppy_basecall_client", shell=True) == 1:
            print("Test stopped. guppy_basecall_client must be in PATH")
            return
        gitDir = os.path.dirname(os.path.realpath(__file__))
        nanoclid = NanoClid(os.path.join(gitDir, "data", "input"), "ADAPTIVE_00", os.path.join(gitDir, "data", "output"), os.path.join(gitDir, "data", "test.bed"), False, "hg19", os.path.join(gitDir,"annotations" ,"gencode.v40lift37.annotation.gtf"), os.path.join(gitDir, "data", "geneList.txt")  , self._SNAKEMAKE_MINION_PARAMETERS, 1, os.environ["HOME"], os.path.join(gitDir, "annotations"), "", os.path.join(gitDir, 'config/config.yaml'),os.path.join(gitDir, 'config/ADAPTIVE.template'), gitDir, os.path.join(os.path.dirname(os.path.realpath(__file__)), 'singularity'), os.path.join(gitDir,"annotations" ,"data"), cpu)
        nanoclid._runNanoClid()
        nanoclid._checkTest(os.path.join(gitDir, "data", "input"), os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'check.template'), os.path.dirname(os.path.realpath(__file__)), cpu)


    def _checkTest(self, inputDir, checkTemplate, gitDir, cpu):
        differences = []
        print("Checking test results ...")
        testDir = os.path.join(gitDir, 'data', 'output')
        if cpu:
            expectedDir = os.path.join(gitDir, 'data', 'expected', 'cpu')
        else:
            expectedDir = os.path.join(gitDir, 'data', 'expected', 'gpu')
        singularityDir = os.path.join(gitDir, 'singularity')
        filesToCheck = open(checkTemplate, "r").readlines()
        filesToCheck = [f.rstrip() for f in filesToCheck]
        for f in filesToCheck:
            if "fastq" in f:
                subprocess.call(f"zcat {f.replace('FOLDER', inputDir)} | sort -Vu > {testDir}/test.fastq", shell=True)
                subprocess.call(f"zcat {f.replace('FOLDER', expectedDir)} | sort -Vu > {expectedDir}/expected.fastq", shell=True)
                diff = subprocess.call(f"diff {testDir}/test.fastq {expectedDir}/expected.fastq", shell=True)
            elif f.split(".")[-1] == "bam":
                diff = subprocess.check_output(f"singularity exec {singularityDir}/bamUtil.sif bam diff --in1 {f.replace('FOLDER', testDir)} --in2 {f.replace('FOLDER', expectedDir)} --all", shell=True)
            elif "pdf" in f:
                if not os.path.exists(f.replace('FOLDER', expectedDir)):
                    print(f"Pipeline failed {f} does not exist")
                    return
            elif os.path.splitext(f)[-1] == ".gz":
                subprocess.call(f"zgrep -v '#' {f.replace('FOLDER', testDir)} > test.tmp", shell = True)
                subprocess.call(f"zgrep -v '#' {f.replace('FOLDER', expectedDir)} > expected.tmp", shell = True)
                diff = subprocess.call(f"diff test.tmp expected.tmp", shell=True)
                subprocess.call("rm test.tmp expected.tmp", shell = True)
            elif os.path.splitext(f)[-1] == ".vcf":
                subprocess.call(f"grep -v '#' {f.replace('FOLDER', testDir)} | cut -f 1,2,3,4,5 > test.tmp", shell = True)
                subprocess.call(f"grep -v '#' {f.replace('FOLDER', expectedDir)} | cut -f 1,2,3,4,5 > expected.tmp", shell = True)
                diff = subprocess.call("diff test.tmp expected.tmp", shell=True)
                subprocess.call("rm test.tmp expected.tmp", shell = True)
            else:
                diff = subprocess.call(f"diff {f.replace('FOLDER', testDir)} {f.replace('FOLDER', expectedDir)} -I '^#'", shell=True)
            if type(diff) == bytes and diff.decode("utf-8") != "":
                differences.append(f.replace('FOLDER', testDir))
            if type(diff) == int and diff == 1:
                differences.append(f.replace('FOLDER', testDir))
        subprocess.call(f"rm {testDir}/test.fastq && rm {expectedDir}/expected.fastq", shell=True)
        print("Checking test results OK")
        if len(differences) > 0:
            print("Some differences have been observed for the following files :")
            for difference in differences:
                print(difference)
            print("Please see the guppy version used for create the expected results on the README.")
            print("Variation can be observed due to guppy version.")
        print("Installation succeed\nNanoCliD is ready to use")

    def _runNanoClid(self):
        if ",".join(self.samples) == "" :
            self.samples = [folder for folder in os.listdir(f"{self.inputFolder}/{self.run}") if self.run in folder and os.path.isdir(f"{self.inputFolder}/{self.run}/{folder}")]
        self.summaryFiles = self.getSpecificPath(f"{self.inputFolder}/{self.run}", "sequencing_summary*.txt", "f", self.samples)
        self.reportFiles = self.getSpecificPath(f"{self.inputFolder}/{self.run}", "report*.md", "f", self.samples)
        self.fast5Dir = self.getSpecificPath(f"{self.inputFolder}/{self.run}", "fast5*", "d", self.samples)
        configFile = self.__createConfig(self.configTemplate, self.outDir, self.run)
        self.runSnakemake(os.path.join(self.gitDir, "workflow/Snakefile"), configFile, self.snakemakeParameters, self.bind, self.cores, self.dryRun)

if __name__ == "__main__":
    parser = ArgumentParser(
	description="Launch NanoClid pipeline")

    subs = parser.add_subparsers(dest = "command")

    install_parser = subs.add_parser("install", help='Install NanoClid images.')
    install_parser.add_argument("-d", "--dir", help = "Folder containing the archives. Please see the README.")

    test_parser = subs.add_parser("test", help="Run unit test")
    test_parser.add_argument("--cpu", action = 'store_true', help = "Run test with a CPU")

    run_parser = subs.add_parser("run", help='Run NanoClid')
    run_parser.add_argument("-b", "--bed", required=True, help="Bed file for adaptive analysis.")
    run_parser.add_argument("-B", "--bind", help="Path to bind. Default is HOME", default=os.environ["HOME"])
    run_parser.add_argument("-c", "--cores", help="Number of cores to be used by snakemake. Default is 1.", default=1)
    run_parser.add_argument("-C", "--containersPath", help="Path to containers tool for NanoClid.", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'singularity'))
    run_parser.add_argument("--cpu", action = 'store_true', help = "Run NanoCliD with CPU")
    run_parser.add_argument("-d", "--dryRun", action='store_true', help="Run nanoclid in dry run mode.")
    run_parser.add_argument("-D", "--snpEffDir", help="Folder containing snpEff annotation files.",  default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'annotations/data'))
    run_parser.add_argument("-f", "--gtf", help="GTF file for circos annotation.", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'annotations/gencode.v40lift37.annotation.gtf'))
    run_parser.add_argument("-g", "--genomeVersion", required=True, help="Genome version for the analysis.")
    run_parser.add_argument("-G", "--gitDir", required=True, help="Path to nanoclid git folder.")
    run_parser.add_argument("-I", "--inputFolder", required=True, help="Path to the folder containing nanopore sequencing data.")
    run_parser.add_argument("-l", "--geneList", required=True, help="File containing genes to look for genomic event in circos representation.")
    run_parser.add_argument("-O", "--outDir", required=True, help="Path to output folder.")
    run_parser.add_argument("-p", "--snakemakeParameters", help="Snakemake parameters to run nanoclid.", default=NanoClid._SNAKEMAKE_MINION_PARAMETERS)
    run_parser.add_argument("-r", "--run", required=True, help="Run ID.")
    run_parser.add_argument("-R", "--refDir", required=True, help="Path to folder containing genome files")
    run_parser.add_argument("-s", "--samples", help="Only run analysis on these samples. Must be comma separated.", default="")
    run_parser.add_argument("-S", "--samplesheet", help="Path to the samplesheet.")
    run_parser.add_argument("-t", "--configTemplate", help="Path to the config template.", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/config.yaml'))
    run_parser.add_argument("-T", "--outputTemplate", help="Path to the output template.", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/ADAPTIVE.template'))

    args = parser.parse_args()

    if args.command == "install":
        NanoClid()._install(args.dir)

    if args.command == "test":
        NanoClid()._runTest(args.cpu)

    if args.command == "run":
        NanoClid(args.inputFolder, args.run, args.outDir, args.bed, args.dryRun, args.genomeVersion, args.gtf, args.geneList, args.snakemakeParameters, args.cores, args.bind, args.refDir, args.samples, args.configTemplate, args.outputTemplate, args.gitDir, args.containersPath, args.snpEffDir, args.cpu)._runNanoClid()
