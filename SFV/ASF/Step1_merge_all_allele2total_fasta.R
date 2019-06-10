#==================================================================================================================================================================
#Step1_merge all sample's allele gene into a total fasta file
#In this script the matrix file named "cgMLST.tsv" and allele be used are required.
#The colnames of the cgMLST.tsv are the final allele be used. Make a txt file named "faname.txt" containing the allele name without ".fasta" 
#Make a floder named using_genes and do: for i in `cat faname.txt`; do cp chewBBACA_step5 output folder/genes_html/${i}_aligned.fasta ./using_genes; done #This is on linux.
#The input data on the R are the cgMLST.tsv and using_genes.
#This script will output a file, suggest to rename the file name, such as total.fasta
#Author: Dongyan Xiong
#==================================================================================================================================================================
#!/usr/bin/Rscript
arg<-commandArgs(T)
setwd(arg[1])
MLST_matrix<-read.table("cgMLST.tsv",header = T,sep = "\t")
MLST_matrix<-as.matrix(MLST_matrix)
rownames(MLST_matrix)<-MLST_matrix[,1]
MLST_matrix<-MLST_matrix[,2:ncol(MLST_matrix)]
dimnames<-list(rownames(MLST_matrix),colnames(MLST_matrix))
MLST_matrix<-matrix(as.numeric(as.matrix(MLST_matrix)),nrow = nrow(MLST_matrix),dimnames = dimnames(MLST_matrix))
library(stringr)
dir1<-dir(arg[2])
setwd(arg[2])
PT3fastaname<-rownames(MLST_matrix)
library(Biostrings)
fa_total<-DNAStringSet()
for (j in 1:ncol(MLST_matrix)){
  fa_col<-DNAStringSet()
  num<-gsub(pattern="[A-Za-z]+",replacement = "",dir1[j])
  num<-gsub(pattern="[0-9]+-",replacement = "",num)
  num<-gsub(pattern="\\.",replacement = "",num)
  Standerdfasta<-readDNAStringSet(dir1[j])
  Matrix_taxno<-MLST_matrix[,j]
  for (i in 1:nrow(MLST_matrix)) {
    faname<-names(Matrix_taxno[i])
    faseqid<-as.numeric(Matrix_taxno[i])
    if(faseqid==0){next;}
    else{
      faseq<-DNAStringSet(toString(Standerdfasta[[faseqid]]))
      names(faseq)<-paste(faname,num,sep = "_")
      fa_col<-c(fa_col,DNAStringSet(faseq))
    }
  }
  fa_total<-c(fa_total,fa_col)
}
out<-tempfile(tmpdir = arg[1])
writeXStringSet(fa_total,out)
