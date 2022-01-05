# RNAInt 2.0
A tool to predict the RNA-interacting residues in a protein sequence.

## Introduction
RNAInt 2.0 is developed for predicting the RNA-interacting residues in the protein using primary sequence information. More information on RNA-interacting residues is available from its web server http://webs.iiitd.edu.in/raghava/rnaint2. This page provide information about standalone version of RNAInt 2.0. Please read/cite RNAInt 2.0 paper for complete information including algorithm behind RNAInt 2.0.<br>

## Standalone
### Minimum Usage
```
python rnaint2.py -i example_seq.fa
```
where example_seq.fa is a input fasta file. This will predict the residues in the seqqeunces provided  in fasta format as RNA-interacting and non-interacting. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

### Complete Usage
```
python dbpred.py [-h] -i INPUT [-o OUTPUT] [-j {1,2,3,4}] [-t THRESHOLD] [-p PATH]
====================================================================================================
Please provide following arguments

  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: Protein sequences in FASTA format
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.6
  -p PATH, --path PATH  Path: Please provide the path of python which has all libraries installed
=====================================================================================================
```
