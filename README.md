# NanoCliD
NanoCliD is a toolkit designed to capture genomic alterations including methylation profile in nanopore data.

NanoCliD uses Snakemake and python.

# Prerequisites

To run NanoCliD you'll need :
  - python3 or higher
  - singularity

# Installation
To install NanoCliD first run :

```python3
python3 nanoclid.py install
```

The script makes each singularity images that will be used by NanoClid and initialises venv to run NanoClid.

Once the installation succeed, you can run the test part to check if NanoCliD is correctly installed. First source virtualenv in venv and run :

```python3
python3 nanoclid.py test -R ${REF_DIR}
```

_REF_DIR must contains fasta file, index files and genome of the genome reference_

If the comparison between run test outputs and expected outputs reports no differences, the installation is completed. 
You can now use NanoCliD.

# Run NanoCliD

To run NanoCliD, the input folder must be organised as follow :


```bash
RUN
|----INJECTION_N
         |----NANOPORE_NOMENCLATURE
                      |barcode_alignment.tsv
                      |duty_time.csv
                      |final_summary.tsv
                      |report.md
                      |sequencing_summary.txt
                      |throughput.csv
                      |----fast5
                      |----other_reports
```

As we used a minion to sequence DNA, this organisation should be generated automatically by the sequencer for each injection.

A minimal call for NanoClid is :

```python3
python3 nanoclid.py run --inputFolder ${INPUT} --run ${RUN} --refDir ${REF_DIR} --genomeVersion hg19|hg38 --bed ${BED} --gitDir ${GIT_DIR}
```
The inputFolder will be set as outputFolder if not set in command line.
Other arguments are available for nanoclid :

- -B, --bind : Path to bind for singularity images. Default is HOME
- -c, --cores : Number of cores used by snakemake (default is 1)
- -C, --containersPath : Path to containers tool for NanoCliD if moved out from the git repository
- -d, --dryRun : Run NanoCliD in dry run mode
- -p, --snakemakeParameters : Snakemake parameters to run NanoCliD (default is "-p --verbose --latency-wait 60 --keep-going --use-singularity")
- -s, --samples : Only run analysis on these samples. Must be comma separated.
- -t, --configTemplate : Path to the config template. (default is config/config.yaml)
- -T, --outputTemplate : Path to the output template. (default is config/ADAPTIVE.template)
