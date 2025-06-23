CHANGES IN VERSION 1.2.0
[ADD]
       o Add GridION
       o Use bin from CNV from BAM for QDNASeq
[BUGFIX]
       o Fix relaunch from pod5/blow5/fast5
       o Deal with empty VCF pepper and empty SV table
       o Get bed even if is in archive folder

CHANGES IN VERSION 1.1.9
[BUGFIX]
       o Fix input folder while looking for input files

CHANGES IN VERSION 1.1.8
[BUGFIX]
       o Fix demultiplexing
       o Harcode genome versions

CHANGES IN VERSION 1.1.7
[BUGFIX]
       o Relaunch from workspace corrected
       o Do not block GUI if MK1C is offline
       o Catch permission denied while cleaning analysis folder
       o QC table with NA value if off target bed is empty
       o SnpEff for hg19 and hg38
       o Remove nanovar

CHANGES IN VERSION 1.1.6
[BUGFIX]
       o Correct sending email if fail  
[ADD]
       o Add an option to set transverse folder

CHANGES IN VERSION 1.1.5
[ADD]
       o Add NanoPlot QC for adaptive and methylation
       o Basecalling can be performed on abacus

CHANGES IN VERSION 1.1.4
[BUGFIX]
       o Deal with fastq from MK1C if demultiplexing already performed
       o Use always newest bed if exists in archive
[ADD]
       o Add .p folder

CHANGES IN VERSION 1.1.3
[BUGFIX]
    o Hotfix : update profile binding, update nanovar

CHANGES IN VERSION 1.1.2
[BUGFIX]
    o Test can be launched on abacus
    o Update path for new genomes
    o Stop launching if renamed run exists in workspace

CHANGES IN VERSION 1.1.1
[BUGFIX]
    o Fix abacus path 

CHANGES IN VERSION 1.1.0
[ADD]
        o Add env/prod environment
        o Add curie/externe environment
        o Add github update facility

CHANGES IN VERSION 1.0.4
[BUGFIX]
        o Handle skip folder from MK1C

CHANGES IN VERSION 1.0.3
[BUGFIX]
        o Calcsub compatibility 

CHANGES IN VERSION 1.0.2
[BUGFIX]
    o Do not check bed nomenclature if methylation run

CHANGES IN VERSION 1.0.1
[BUGFIX]
    o Remove bamstat file
    o Correct transverse error mail
    o Methylation file format for FMP export
    o Correct bed if nomenclature issue
NEW FEATURES
    o Add gitlab templates

CHANGES IN VERSION 1.0.0
NEW FEATURES
    o Init pipeline
