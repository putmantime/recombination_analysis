# recombination_analysis

Python Class for determinating recombination break points of child recombinant by comparing to parents via a sliding window blast approach

Usage:
#1 need to create a blastdb from multi.fasta (unaligned) of parents

e.g.
parents.fasta

>parent1
gcttgcgtgcgcgtcgtcgtc....
>parent2

makeblastdb -in parentes.fasta -parse_seqids -dbtype nucl

#2 rec_pipeline.py
Usage: provide 6 command line arguments

$ rec_pipeline.py child_genome window_size (integer 1000) step_size (integer 100) blastdatabase parent1_name (name of fasta header in blastdb file without ">") parent2_name (name of fasta header in blastdb file without ">")

example:  
./rec_pipeline.py sequences/rec_genome_test.fasta 1000 100 sequences/multi_parent_test.fasta parent_test_1 parent_test_2


Output: 
Table of genome coord, window name, parent code (int)


|coord| window| score|
| ------------- | ------------- | ------------- |  
|1980| rec_child_test_window_2_coord_1980 |0|
|5940| rec_child_test_window_6_coord_5940 |0|
|2970| rec_child_test_window_3_coord_2970 |0|
|990 |rec_child_test_window_1_coord_990   |2|
|4950| rec_child_test_window_5_coord_4950 |0|
|6930| rec_child_test_window_7_coord_6930 |0|
|9900| rec_child_test_window_10_coord_9900| 1|
