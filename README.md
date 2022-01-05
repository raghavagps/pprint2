# RNAInt 2.0
A tool to predict the RNA-interacting residues in a protein sequence.

## Introduction
RNAInt 2.0 is developed for predicting the RNA-interacting residues in the protein using primary sequence information. More information on RNA-interacting residues is available from its web server http://webs.iiitd.edu.in/raghava/rnaint2. This page provide information about standalone version of RNAInt 2.0. Please read/cite RNAInt 2.0 paper for complete information including algorithm behind RNAInt 2.0.<br>

## Standalone

### Important Instruction
In order to run the code, you need to have model file which is stored as prog.zip. Kindly unzip it using the following command:
```
unzip prog.zip
```

### To get the help
```
python rnaint2.py -h
```
### Minimum Usage
```
python rnaint2.py -i example_seq.fa
```
where example_seq.fa is a input fasta file. This will predict the residues in the seqqeunces provided  in fasta format as RNA-interacting and non-interacting. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

### Complete Usage
```
python dbpred.py [-h] -i INPUT [-o OUTPUT] [-j {1,2,3,4}] [-t THRESHOLD] [-p PATH]
```
```
====================================================================================================
Description of each argument:

  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: Protein sequences in FASTA format
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.6
  -p PATH, --path PATH  
                        Path: Please provide the path of python which has all libraries installed
=====================================================================================================
```
* Input File: It allow users to provide input in i) FASTA format (standard) and ii) Simple Format (Single line). The program will generate the patterns of length 17 for each sequence.
* Output File: Program will save result in CSV format, in case user do not provide output file name, it will be stored in outfile.csv. In the output file, three lines will be displayed for each sequence, where, first line represents the ID, second row explains the sequence, and third line exhibits the '+' and '-' sign, where '+' corresponds to residues which are involved in RNA-interaction, and '-' represents the residues which are not involved in RNA-interaction.
* Threshold: User should provide threshold between 0 and 1, please note score is propotional to the RNA interacting potential of the resiudes.



