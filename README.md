# NanoCliD
NanoCliD is a toolkit designed to capture genomic alterations including methylation profile in nanopore data.

# Prerequisites

To run NanoCliD you'll need :
  - python3 or higher
  - singularity

# Installation
To install NanoCliD first run :

python3 nanoclid.py install

Once the installation succeed, you can run the test part to check if NanoCliD is correctly installed :

python3 nanoclid.py test -R ${REF_DIR}

*REF_DIR must contains fasta file, index files and genome of the genome reference
