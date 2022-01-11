#!/usr/bin/bash

#Trim the low-quality data
fastp -i data_R1.fq.gz -I data_R2.fq.gz -o data_clean_R1.fq.gz -O data_clean_R2.fq.gz

#Assemble the genome
megahit -t 25 -m 0.8 -1 data_clean_R1.fq.gz -2 data_clean_R2.fq.gz -o Assembled-genome-data

#Full genome construction
#1.Add the total of contigs into the windows based software SeqMan to archive the longer scaffolds.
#2.Using the BLASTN to identify the A.baumanii sequences. And full the gap with PCR-sequencing.
