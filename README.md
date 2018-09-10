# qiaseq-rna
```
usage: python3 qiaseq-rna/run.py --runList runList.tsv [--numCores] [--ge] [--fc] [--vc]

speRna pipeline: perform primer QC, gene expression quantification, fusion detection, and SNP/INDEL calling

optional arguments:

    -h, --help            show this help message and exit
    --numCores            number of threads to be used, default = 1
    --ge                  perform gene expression quantification
    --fc                  perform fusion calling
    --vc                  perform snp/indel calling

The input table, e.g. runList.tsv, is a tab-delimited file with at least 10 columns:

    0   runpath
    1   runname
    2   samplename
    3   speSide
    4   r1universal
    5   r2universal
    6   umiLen
    7   umiOffset
    8   panelId
    9  platform
    
Lines starting with "#" and columns beyond the 10th are ignored. The pipeline expects fastq or fastq.gz files as:

    /runpath/runname/samplename*_R1*.fastq*
    /runpath/runname/samplename*_R2*.fastq*

for Illumina reads and only the first one for Ion reads. The platform is defined as:

    0	illumina (except nextseq)
    1	nextseq
    2   ion torrent  
```


