# Seq2Karyotype

Seq2Karyotype is a tool designed for copy number alterations analysis and karyotyping using bulk whole-genome sequencing data. This README provides instructions for installation, configuration, and usage of Seq2Karyotype, along with examples and helpful tips.

## Change Log

## Installation

Download the package, then follow the instructions below.

```
# setup the virtual environment
conda -create -n s2k python=3.11
# activate the virtual envrionment
conda activate s2k
# install the package
cd /path/to/download/folder
pip install Seq2Karyotype
```

## Download Auxiliary Files

CytoBand annotation are provided in the folder `aux`. You will have to download the SNP whitelist yourself. Both lists for hg19, hg38 are provided under https://drive.google.com/drive/folders/1ontR5R_dfl8uJc50PwHA2hxEBIcdkmn7?usp=drive_link. You can use your own SNP list if you wish, just make sure to follow the same format as in the provided files.

## Usage 

1. Convert your variant calling file to allele count format.

An example snippet is provided below

```
#!/bin/bash

filename=$( echo $1 | cut -d'/' -f12 | cut -d'.' -f 1 )
echo $filename

zcat -1 $1 | grep -v -e '#' | grep AD | awk 'BEGIN {print "Chr\tPos\tRef\tAlt\tType"};{split($10,a,":"); split (a[2],b,","); if (length(b) == 3) print $1"\t"$2"\t"b[1]"\t"b[2]"\tSNP"}' > $filename.ac
```

The output allele count file should contain at least five columns

| chromosome | position | reference_count | alternative_count | type |
|----------|------|----------|-------|-------|
chr1|10513|12|3|SNP
chr1|13868|1|6|SNP
chr1|14464|7|35|SNP

2. Get a local copy of the configuration file.

```
s2k getconfig
```

3. Edit the `[Input]` and `[InputColumns]` sections in the configuration file to suit your data. We don't recommend changing additional parameters when starting out.

4. Run S2K to generate output files.

```
s2k analyze -i sample_name.txt -c your_configuration_file.ini -s sample_name -m0 diploid_coverage -mc merging_coeff
... Lots of log info! ...
12:18:19 S2K.WGS: INFO: Ready to report!
All done
```

Note that `m0` and `mc` are optional. `m0` will be calcuate by the diploid detection algorithm, and can be adjusted if incorrect. `mc` is the merging coefficient that controls the merge between segements.

There will be 4 output files in total:

`sample_name.bed`: The main output file contains the segmentation information.

`sample_name.dat.gz`: The output for individual chromosome information.

`sample_name.log`: Your good old log file.

`sample_name.par`: The parameters used including reference diploid coverage, karyotype models, etc,.

5. Launch the viewer for visualization.

```
# if you are running this on your local machine, such as a desktop or laptop:
s2k viewer

# if you are running this on a remote machine, like a cluster:
s2k viewer --remote
**********
Access dashboard in browser via: http://10.220.16.129:39738
**********
INFO:     Started server process [11716]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://10.220.16.129:39738 (Press CTRL+C to quit)
```

Load the output files

![example](./images/example_file_selection.gif)

## Acknowledgments

If you find this tool useful, please star this repo and cite our paper (when it comes out) :)

Limeng Pu et al. "Seq2Karyotype (S2K): A method for deconvoluting heterogeneity of copy number alterations using single-sample whole-genome sequencing data" [abstract]. Proceedings of the American Association for Cancer Research Annual Meeting 2024; Part 1 (Regular Abstracts); 2024 Apr 5-10; San Diego, CA. Philadelphia (PA): AACR; Cancer Res 2024;84(6_Suppl):Abstract nr 7419

This README file is written by Limeng Pu.


