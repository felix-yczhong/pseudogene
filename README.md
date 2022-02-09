# Pseudogene Copy Number Analysis Tool
## Table of Content
[Introduction](#Introduction)
[Installation](#Installation)
[Configuration](#Configuration)
[Usage](#Usage)

## Introduction
pseudogene copy number analysis tool compares corresponding bases of a true gene and its pseudo gene(s) to calculate their copy numbers.

## Installation
Use the following commands to download source code, create and activate virtual environment, and install the program

```
$ git clone
$ cd pseudogene
$ python3 -m venv pseudogene_venv
$ source pseudogene_venv/bin/activate
$ pip3 install --editable=.
```

## Configuration
Modify config.json to change program configuration. Function of each field is explained below.

Note that all gene sequence coordinate must be 0-based.

genes: a list of gene group on which this program executes
output: output file name format

* temp_files: file name format used during exeuction, only used by developer
    * vcf_path: file name format used to create a vcf file
    * nirvana_path: file name format used for Nirvana to create its summary

* tools: path of tools
    * fasta_loc: path of reference genome fasta, used to query reference sequence
    * nirvana_path: path of Nirvana program

* seq_align_map: file name formats of sequence alignment maps of each gene, where the program outputs only when their corresponding reference base differs
    * seq_align_map files must be save in seq_align_maps/

* force: file name formats of sequence alignment maps of each gene that forces the program to output
    * force files must be saved in seq_align_maps/

* gene group:
    * true gene: name of true gene
    * pseudo genes: a list of pseudo genes
    * NM_number: NM number used to query amino acid change
    * true gene region: gene coordinates of true gene
    * pseudo gene regions: gene coordinates of pseudo genes
    * CNV ratio range: a 2-tuple, e.g. (I, J)
        if I and J are integers, then only search ratio N/M where N in {0, ... , I} and M in {0, ... , J}
        if I and J are lists of integers, then only search ratio N/M where N in I and M in J
        I and J can be either a list or an integer.

## Usage
To get a list of arugments and options.
```
$ pseudogene --help
```
Arguments:
  BAM_LIST: a list of bam file paths

Options:
  --ncpus INTEGER              number of cores to use  [default: 1]
  -r, --reference [hg19|hg38]  reference genome that was used for alignment
                               [default: hg38]
  --debug                      write out scale factors and all sequence alignment details instead of just summary. Do NOT run debug mode with large number of genes.
  -o, --output DIRECTORY       specify output directory.
  --help                       Show this message and exit.

This program uses min(ncpus, number of genes) processes. Each process uses ncpus threads.

## Intrepretation
