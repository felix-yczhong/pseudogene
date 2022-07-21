# True Gene and Pseudogene Copy Number Analysis Tool
## Table of Content
[Introduction](#introduction)

[Installation](#installation)

[Usage](#usage)

[Output](#output)

[Configuration](#configuration)

[Clinical Interpretation](#clincal-interpretation)

## Introduction
True Gene and **P**seudogene **C**opy **N**umber **A**nalysis **T**ool (PCNAT) compares differentiating bases at aligned positions betwen a target true gene and pseudogene(s) to estimate their copy numbers.

PCNAT contains 3 main functions:
1. build_db: build a databse of genomic information from .gff3 file
2. query_alignment: align a target true and pseudogene(s) and find positions where their reference bases differ
3. run: run the algorithm to estimate their copy numbers

## Installation
Use the following commands to download source code, create and activate virtual environment, and install the program

```
$ git clone
$ cd pseudogene
$ python3 -m venv pseudogene_venv
$ source pseudogene_venv/bin/activate
$ pip3 install --editable=.
```

## Usage
To get a list of arugments and options.
```
$ pseudogene --help
Usage: pseudogene [OPTIONS] COMMAND [ARGS]...
Options:
    --version  Show the version and exit.
    --help     Show this message and exit.

Commands:
  build_db         Build database for sequence alignment
  query_alignment  Query sequence alignment.
  run              Run Pseudogene calculation tool.

```
```
$ pseudogene build_db --help
Usage: pseudogene build_db [OPTIONS] DATABASE_SOURCE

  Build database for sequence alignment

Options:
  -r, --reference [hg19|hg38]  [default: hg38]
  --database_path PATH         path where the database will be built. Will use
                               the one in config.json if empty.
  --help                       Show this message and exit.
```
```
$ pseudogene query_alignment --help
Usage: pseudogene query_alignment [OPTIONS]

  Query sequence alignment.

Options:
  --debug                      output aligned sequence in sam format
  -r, --reference [hg19|hg38]  [default: hg38]
  --gene_groups TEXT           query only for specific gene group(s), gene
                               group name(s) must be in config.json
  --help                       Show this message and exit.
```
```
$ pseudogene run --help
Usage: pseudogene run [OPTIONS] BAM_LIST...

  Run Pseudogene calculation tool.

Options:
  --ncpus INTEGER              number of cores to use.  [default: 16]
  -r, --reference [hg19|hg38]  reference genome used for alignment  [default:
                               hg38]
  --debug                      write out details such as scale factors,
                               control genes average coverage and true gene &
                               pseudogene average coverage.
  --profile_path FILE          specify output path of profiling result, not
                               used if profiling is off  [default:
                               current directory]
  --profile                    enable profiling.
  -o, --output DIRECTORY       specify output directory.  [default:
                               current directory]
  -c, --control FILE           takes a file path which contains all control
                               samples' file paths
  --help                       Show this message and exit.
```

This program uses min(ncpus, number of genes) main processes for each gene. Each main process uses ncpus sub-processes to handle input bam files and estimate copy numbers.

## Output
Each case has one output file. If debugging mode is turned on, control genes average coverage, target true gene and pseudogene(s) average coverage and scale factors for all cases will be output as well.
The meaning of each output field is explained below
* ratio
    * exception
        * is_exception: whether this position is specified as exception pos in config.json
    * ref: reference bases for each position
    * bases: read counts for reference bases for each positoin
    * copy number: estimated copy number for true gene and pseudogene(s)
* sum
    * bases: sum of each base over true gene and pseudogene(s) at each aligned position.
* true: true gene's genomic information such as exon, chromosome, position, reference base, read counts of each base
* pseudo: each pseudogene's genomic information such as exon, chromosome, position, reference base, read counts of each base

## Configuration
Modify config.json to change program configuration. Function of each field is explained below.

### Main Configuration
* gene_group: the gene group(s) you wish to run on, by default the gene group name is the name of its assocaited disease
* output: the name of the output files. PCNAT only outputs in XLSX format
* debug_output: the name of the debugging information file, which contain information such as each control genes' average coverage, target true gene and pseudogene(s) average coverage, etc. PCNAT only outputs in .tsv format.
* summary_tab: tab names of the output files

* tools: location of tools used by PCNAT
    * fasta_loc: location of reference genome fasta, used to query reference sequence
    * BLAST+_path: location of BLAST+ program, used to find alignment relations between a target true gene and pseudogene(s)
    * BLAST+_command: command used to call BLAST+ program, by default uses megablast for highly similar sequences
    * nirvana_path: location of Nirvana program, used to find associated amino acid changes

* data: location of gene sequences, exon sequences and alignment relations
    * database_path: location of genomic information database
    * fasta:
        * true_seq: location of a true gene's sequence
        * pseudo_seq: location of a pseudogene's sequence
        * true_exon: location of a true gene's exons
        * all_true_exon: location of all true genes' exons
        * all_pseudo_seq: location of all pseudogene(s)' sequences
    * sam:
        * aligned_seq: location of alignment relation between a target true gene and pseudogene(s)

To manually view the alignment relations, open IGV, load all_true_exon as reference genome and open aligned_seq file.

### Gene Group Configuration
All genome coordinate must be 0-based. For each gene group, in general the format is demonstrated as the following.

Please also see default config.json for examples.

true_gene should contain exactly 1 string. More than 1 true gene is not supported.

pseudo_gene can contain more than 1 strings. More than 1 pseudogene is supported.
```
"gene_group_name_1": 
{
    "true_gene": ["true_gene_name_1"],
    "pseudo_gene": ["pseudogene_name_1"],
    "hg38": 
    {
        "true_gene_name_1": {
            "NM_number": "true_gene_name_1_NM_number",
            "gene_id": "true_gene_name_1_gene_id",
            "transcript_id": "true_gene_name_1_transcript_id",
            "exception": null or list of list of 3 elements
        },

        "pseudogene_name_1": {
            "gene_id": "pseudogene_name_1_NM_number",
            "transcript_id": "pseudogene_name_1_transcript_id",
            "exception": null or list of list of 3 elements
        },
    }      
}
```
Each field is explained below.

| Name | Data Type | Usage |
|:---:|:---:|:---:|
| gene group | string | name of the gene group to query, used in query_alignment function |
| true_gene | list of string | name of the true gene, or gene used as reference template in alignment |
| pseudo_gene | list of string | name of the pseudogene, or other genes |
| NM_number | string | accession number for any gene, used to look up amino acid change by Nirvana |
| gene_id | string | accession number for any gene, used to look up gene genomic range in GENCODE database |
| transcript_id | string | accession number for any gene, used to look up gene genomic range in GENCODE database |
| exception | Null or list of lists of string  | Used to manually call out positions that have same reference base. The inner list must contain exactly 3 elements, which are chromosome, position, and alternate base. Alternate base is used to indicate to the program which base should be used to calculate copy number for pseudogene(s) instead of using reference base. |

## Clincal Interpretation
For SMA, cases with fewer than 2 copies of SMN1 are considered patients; cases with 2 copies of SMN1 but fewer than 2 copies of SMN2 are considered silent carriers.

For Gaucher disease, cases with fewer than 2 copies of GBA are considered patients; cases with 2 copies of GBA but fewer than 2 copies of GBAP1 are considered silent carriers.
