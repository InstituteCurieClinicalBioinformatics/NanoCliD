# NanoCliD
The Clinical Bioinformatics Team at Curie Institute in Paris provides NanoCliD (Nanopore Clinical Diagnosis), a toolkit designed to capture genomic alterations including methylation profile in nanopore adaptive sampling data.

NanoCliD uses Snakemake and python.

# Prerequisites

To run NanoCliD you'll need :
  - python > 3.7 and pip3
  - python3-venv
  - singularity
  - [guppy](https://community.nanoporetech.com/downloads)

# Installation
To install NanoCliD first run :

```python3
python3 nanoclid.py install
```

The script will download from google drive singularity archives, input and expected test files and annotations archives for hg19 that will be used by NanoClid and initialises venv to run NanoClid.

However, due to high size of the archives and google politics, the download can fail returning "Too many requests" error. In this case, just wait and try later.
You can also manually download the archives at the following links :

singularity images part1 : https://drive.google.com/file/d/1rnOcr8M_lH3KSu3dIKpi6C06cxSt3zMN/view?usp=share_link

singularity images part2 : https://drive.google.com/file/d/1-ZkxulbbVv54y__Zz6y-Xahh-YwftPUm/view?usp=share_link

annotations : https://drive.google.com/file/d/1RQa5QHgMDmcFR-IRSKwbUiMn_QMMoYKc/view?usp=share_link

input test files : https://drive.google.com/file/d/1LpguQ1aUFQPGcnh4iow73rExaHT-E96_/view?usp=share_link

expected test files : https://drive.google.com/file/d/1hkn4YxjmID30Rqil58XjdiEUzm_1RTLq/view?usp=share_link


If you downloaded manually the archives, please launch the following command to install :

```python3
python3 nanoclid.py install -d ${DOWNLOAD_DIR}
```
Where DOWNLOAD_DIR is the folder containing the archives which must be different from the git folder.

Once the installation succeed, you can run the test part to check if NanoCliD is correctly installed. First source activate in venv and run :

```bash
source venv/bin/activate
python3 nanoclid.py test
```

/!\ : if you do not have a GPU, you can run the test using CPU with the --cpu option.
However, we want to remind you that adaptive sampling needs a GPU. You will be able to run the test part with a CPU but not the run part.

The CPU expected results have been obtained using guppy 6.3.8+d9e0f64.
The GPU expected results have been obtained using guppy 5.0.17+99baa5b27.
The comparison between test run outputs and expected outputs could report differences if your guppy version differs from the reported ones.
If the test runs without error, the installation is completed.
You can now use NanoCliD.

# Run NanoCliD

To run NanoCliD, the input folder must be organised as follow :


```bash
RUN
|----INJECTION_N
         |----NANOPORE_NOMENCLATURE
                      |report.md
                      |sequencing_summary.txt
                      |----fast5

# for the small input test dataset
ADAPTIVE_00
├── ADAPTIVE_00_1
│   └── subdir1
│       └── subdir2
│           ├── fast5
│           │   ├── FAK42335_2bf4f211a2e2d04662e50f27448cfd99dafbd7ee_0.fast5
│           │   └── FAK42335_2bf4f211a2e2d04662e50f27448cfd99dafbd7ee_100.fast5
│           ├── report.md
│           └── sequencing_summary.txt
└── ADAPTIVE_00_2
    └── subdir1
        └── subdir2
            ├── fast5
            │   └── FAK42335_2bf4f211a2e2d04662e50f27448cfd99dafbd7ee_200.fast5
            ├── report.md
            └── sequencing_summary.txt

```

As we used a minion to sequence DNA, this organisation should be generated automatically by the sequencer for each injection.

A minimal call for NanoClid is :

```python3
python3 nanoclid.py run --inputFolder ${INPUT} --run ${RUN} --refDir ${REF_DIR} --genomeVersion hg19|hg38 --bed ${BED} --gitDir ${GIT_DIR}
```
The inputFolder will be set as outputFolder if not set in command line.
$REF_DIR must contains fasta file, index files and genome file of the reference genome.

Other arguments are available for nanoclid :

- -B, --bind : Path to bind for singularity images. Default is HOME
- -c, --cores : Number of cores used by snakemake (default is 1)
- -C, --containersPath : Path to containers tool for NanoCliD if moved out from the git repository
- -d, --dryRun : Run NanoCliD in dry run mode
- -p, --snakemakeParameters : Snakemake parameters to run NanoCliD (default is "-p --verbose --latency-wait 60 --keep-going --use-singularity")
- -s, --samples : Only run analysis on these samples. Must be comma separated.
- -t, --configTemplate : Path to the config template. (default is config/config.yaml)
- -T, --outputTemplate : Path to the output template. (default is config/ADAPTIVE.template)

Results tree of NanoCliD

```bash

├── ADAPTIVE_00_1
│   ├── FASTQ
│   │   └── ADAPTIVE_00_1.fastq.gz
│   ├── Mapping
│   │   ├── ADAPTIVE_00_1.bam
│   │   └── ADAPTIVE_00_1.bam.bai
│   ├── Methylation
│   │   ├── ADAPTIVE_00_1.index.txt
│   │   └── ADAPTIVE_00_1.methylation_calls.tsv
│   └── QC
│       ├── ADAPTIVE_00_1_mergedStats_offtarget.txt
│       └── ADAPTIVE_00_1_mergedStats_ontarget.txt
├── ADAPTIVE_00_2
│   ├── FASTQ
│   │   └── ADAPTIVE_00_2.fastq.gz
│   ├── Mapping
│   │   ├── ADAPTIVE_00_2.bam
│   │   └── ADAPTIVE_00_2.bam.bai
│   ├── Methylation
│   │   ├── ADAPTIVE_00_2.index.txt
│   │   └── ADAPTIVE_00_2.methylation_calls.tsv
│   └── QC
│      ├── ADAPTIVE_00_2_mergedStats_offtarget.txt
│      └── ADAPTIVE_00_2_mergedStats_ontarget.txt
└── All
    ├── Circos
    │   ├── ADAPTIVE_00_gene1.csv
    │   ├── ADAPTIVE_00_gene2.csv
    │   └── circos.txt
    ├── CNV
    │   ├── ADAPTIVE_00_amp_plot.pdf
    │   ├── ADAPTIVE_00_cnv_plot.pdf
    │   ├── ADAPTIVE_00_cnv_plot_segments.txt
    │   ├── ADAPTIVE_00_gc.pdf
    │   ├── ADAPTIVE_00_readlength.pdf
    │   └── ADAPTIVE_00_readlength_summary.txt
    ├── Mapping
    │   ├── ADAPTIVE_00.bam
    │   └── ADAPTIVE_00.bam.bai
    ├── Methylation
    │   └── ADAPTIVE_00.methylation_calls.tsv
    ├── QC
    │   ├── ADAPTIVE_00_offTarget.bed
    │   ├── ADAPTIVE_00_qc_offtarget.tsv
    │   ├── ADAPTIVE_00_qc_ontarget.tsv
    │   ├── ADAPTIVE_00_run_qc_offtarget.tsv
    │   └── ADAPTIVE_00_run_qc_ontarget.tsv
    ├── SNV
    │   ├── ADAPTIVE_00_annot.vcf
    │   ├── ADAPTIVE_00.clair3.vcf.gz
    │   ├── ADAPTIVE_00.clair3.vcf.gz.tbi
    │   ├── ADAPTIVE_00.nanocaller.vcf.gz
    │   ├── ADAPTIVE_00.nanocaller.vcf.gz.tbi
    │   ├── ADAPTIVE_00.pepper.vcf.gz
    │   ├── ADAPTIVE_00.pepper.vcf.gz.tbi
    │   ├── ADAPTIVE_00.vcf.gz
    │   └── ADAPTIVE_00.vcf.gz.tbi
    └── SV
        ├── ADAPTIVE_00.annotSV.tsv
        ├── ADAPTIVE_00.merged.annotated.tsv
        ├── ADAPTIVE_00.merged.vcf
        ├── cuteSV
        │   └── ADAPTIVE_00.cuteSV.vcf
        ├── Nanovar
        │   └── ADAPTIVE_00.nanovar.pass.vcf
        ├── sniffles
        │   └── ADAPTIVE_00.sniffles.vcf
        └── SVIM
            └── ADAPTIVE_00.svim.vcf

```

In case of troubleshooting please contact : bioinfo-clinique[AT]curie.fr



# References
[JOBIM 2022](https://jobim2022.sciencesconf.org/resource/page/id/18): **Improving clinical diagnosis using Nanopore Adaptive Sampling and NanoCliD.**
Eleonore Frouin, Kevin Merchadou, Mathilde Filser, Abderaouf Hamza, Elodie Girard, Nicolas Servant, Julien Masliah-Planchon and Victor Renault
