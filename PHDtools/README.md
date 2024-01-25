# PHDtools: Pathogen High-dimensional Data tools for rapid identify the potential pathogen and multi-dimensionally analyze the genetic variations
![](https://img.shields.io/badge/Platform-Linux|Ubuntu(16.04~22.04)-green)
![](https://img.shields.io/badge/Install_with-|github-orange)

**The PHDtools platform is a graphical user interface (GUI) based intranet online platform. Computers in the same LAN could access this platfrom by visiting the intranet IP address of the work station which installed with PHDtools.**

```
  ____  _   _ ____  _              _     
 |  _ \| | | |  _ \| |_ ___   ___ | |___ 
 | |_) | |_| | | | | __/ _ \ / _ \| / __|
 |  __/|  _  | |_| | || (_) | (_) | \__ \
 |_|   |_| |_|____/ \__\___/ \___/|_|___/
```                                      

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Required dependencies](#required-dependencies)
- [Usage](#usage)
- [License](#license)
- [Citation](#citation)

## Introduction
**PHDtools is an intranet online platform that can process the mNGS/NGS data of pathogen genomes. This platform is helpful for the labs and the frontliner workers who have limited background in bioinformatics. Briefly, PHDtools have a total of 15 functions described as follows.**

Core functions:
```
1. Pathogen X Identification (mNGS data analysis)
2. Target Pathogen Genome Assembly (De novo assemble the pathogen genome)
3. Intra-host Variation Identification (Identify within-host minor alleles)
4. Inter-host Variation Identification (Identify mutations including SNP and Indel between two isolates)
5. Lineage-level Variation Identification (Identify common mutations at the lineage-level of the pathogen)
```
Other functions:
```
1. Reverse or Complement of the Sequences
2. Translate Nucleotide Sequences to Amino Acide Sequences
3. Extact the Target Sequences
4. Sequences Alignment
5. Genome Annotation
6. Evaluation of Sequencing Quality
7. Phylegenetic Strain Typing
8. Genome Annotation
9. Identification of Conserved Regions of the Pathogen
10. PCR-primers Mapping
```

## Installation
**The PHDtools platform was installed and tested successfully on Ubuntu 16.04/18.04/20.04/22.04 LTS systems. This platform relies on the LAMP (Linux/Apache/MySQL/PHP/Perl) based Web sever. The essential softwares are recommended to be installed by the follows description.**

- __Kindly remind: if you have any trouble in the installation, kindly please contact hpwei@wh.iov.cn, yujp@wh.iov.cn and xiongdongyan18@mails.ucas.ac.cn, we will arrange some online assistance for you to complete the installation of PHDtools__

- __First: Install the Apache software__

```bash
#install apache2
$sudo apt-get install apache2

#start apache2 server
$sudo /etc/init.d/apache2 start
```

Check Apache server by visiting http://127.0.0.1/index.html

- __Second: Install the MySQL software__

```bash
#install MySQL
$sudo apt-get install mysql-server
$sudo apt-get install mysql-client
$sudo apt-get install libmysqlclient-dev
```

- __Third: Install the PHP software__

```bash
$sudo apt-get install software-properties-common
$sudo add-apt-repository ppa:ondrej/php && sudo apt-get update
$sudo apt-get -y install php7.2
```

Check PHP by visiting http://127.0.0.1/index.php

- __Fourth: Install the PHDtools__

```bash
#download all codes from current github page
$sudo mv -f html/ /var/www/
$sudo mv -f other_scripts/ /var/www/
#mkdir a set of folders as follows to establish the essential paths of PHDtools in your ~/ or home folders
Web-Server
├── AlignPrimers
├── Basic-operation
│   ├── align-seq
│   ├── extract-seq
│   ├── reverse-complen
│   └── translate-seq
├── Compare-genome
├── ConserveSeq
├── Genome-annotation
├── Genome-assemble
├── Microbiome
├── MinorAllele
├── Primer-design
├── Result-store
│   ├── AlignPrimers
│   ├── Compare-genome
│   ├── ConserveSeq
│   ├── Genome-annotation
│   ├── Genome-assemble
│   ├── Microbiome
│   ├── MinorAllele
│   ├── Primer-design
│   ├── Strain-typing
│   ├── Uniq-mut
│   └── WGS-quality
├── Strain-typing
├── Uniq-mut
└── WGS-quality
```

Edit the path within the scripts from /var/www/html/ and /var/www/other_scripts/
Edit the path to required softwares within the scripts /var/www/html/ and /var/www/other_scripts/

## Required dependencies
The dependencies are recommended to be installed with __[Bioconda](https://bioconda.github.io/index.html)__
- [Perl](http://www.perl.org/get.html)
- [perl-bioperl](http://metacpan.org/pod/BioPerl)
- [fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqscreen)
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [NCBI-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [megahit](https://github.com/voutcn/megahit)
- [samtools](https://github.com/samtools/samtools)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Prokka](https://github.com/tseemann/prokka)
- [prodigal](https://github.com/hyattpd/Prodigal)

## Usage
**All the computers in the same LAN can visit the PHDtools platform via intranet IP address of the work station installed with PHDtools. Click into the page of corresponding function and upload the raw data to PHDtools. The results will be release later.**

## License
**Academic use is free of charge. For commerical use, please contact hpwei@wh.iov.cn, yujp@wh.iov.cn and xiongdongyan18@mails.ucas.ac.cn for permission.**

## Citation
**Our article is in preparation, https://github.com/MisgaXiong/Sripts_for_Paper/edit/master/PHDtools can be cited before our paper is published。**
