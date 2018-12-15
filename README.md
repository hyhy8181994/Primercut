# Primercut 
*version 0.1*

A simple python script to cut the primers (forward and reverse) from input fasta file. The input primer can have ambiguous base (M and N). If the program does not detect the primers in the first 50 bases and last 50 bases of sequence with selected matching criteria, it will discharge the sequence to unk file. 


## Options

__-in_file__  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the name of the input FASTA file

__-out_file__ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the name of the trimmed reads file

__-unk_file__ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the name of the file containing unprocessed reads

__-n_mismatch__ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the tolerance for mismatches in the forward or reverse primer sequence

__-min_len__ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the minimum length a sequence must be in order for it to be processed

__-forward__ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the forward primer

__-reverse__ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the reverse primer

## Usage Example 

python3 &nbsp; primercut.py &nbsp; -in_file paired.fasta &nbsp; -out_file out.fasta &nbsp; -unk_file unk.fasta &nbsp; -n_mismatch 3 -min_len 150 &nbsp; -forward GTGCCAGCMGCCGCGGTAA &nbsp; -reverse ACAGCCATGCANCACCT
