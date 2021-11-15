# A tutorial for identification of unique mutations of SARS-CoV-2 variants
In order to establish a simple and accurate PCR based approach for determining the SARS-CoV-2 VOCs, we developed a pipeline to determine the unique mutations for SARS-CoV-2 variants.

Author: Dongyan Xiong
Prof. Hongping Wei's Lab

Here is the usage of the pipeline. All scripts are recommended to run on the Linux system.

### Usage
1.Download the metadata information table of each SARS-CoV-2 genome from the 2019 novel coronavirus resource database (https://ngdc.cncb.ac.cn/ncov/release_genome). Then, rename the downloaded file as metadata.tsv;

2.Construct a high-quality SARS-CoV-2 genome information table using the command line:
```
perl HighHumanStrain.pl metadata.tsv > GISAID_HighHomoStrain_metadata.tsv
```

3.Extract the accession ids from the constructed high-quality SARS-CoV-2 genomes information table using the command line: 
```
perl splitGISAID_ACC.pl GISAID_HighHomoStrain_metadata.tsv
```

4.Batch download the high-quality SARS-CoV-2 genomes from GISAID database.

5.Merge all downloaded genome files (.fasta format) using the command line:
```
perl mergeseq.pl
```

6.Rename the merged genome file; for example, GISAID_HighQC_Genomes.fasta.

7.Merge the genomes and the corresponding information using the command line:
```
perl strainlineage_stat.pl GISAID_HighHomoStrain_metadata.tsv GISAID_HighQC_Genomes.fasta &> run.log
```
The command will output a file named Match_xxx_Lineage.tsv

8.With the output files “Match_xxx_Lineage.tsv” and “GISAID_HighQC_Genomes.fasta”, the core analysis can be run using the command line:
```
perl batchlin.pl Match_xxx_Lineage.tsv and GISAID_HighQC_Genomes.fasta &> run.log
```

9.To obtain the mutation profile within each lineage, run the command line:
```
perl lineageVarWhole.pl Total_strainVar.tsv 0.5 100
```
The parameters 100 and 0.5 means: (1) There must be at least 100 genomes in one lineage. (2) Only mutations present in over 50 % of the genomes within a lineage are considered as the common mutations in this lineage.

10.To obtain the mutation relationship among lineages, run the command line:
```
perl lineageWholeVarstat.pl sig_LineageVar.tsv > results.tsv
```
