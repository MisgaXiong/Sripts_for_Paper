#!/usr/bin/perl -w

=head1 Description

	Comparative sequences analysis
	
Usage

	Genome vs Genome (without amino acid annotation):
		perl compare-genome.pl -WorkPath <path> -datatype <fragment> -seqtype <genome> -reference <reference genome> -isolate <isolate genome> -StorePath <output path>
	
	Genome vs Genome (with amino acid annotation):
		perl compare-genome.pl -WorkPath <path> -datatype <fragment> -seqtype <genome> -reference <reference genome> -isolate <isolate genome> -gff <gff file> -codon <codon table> -StorePath <output path>
		
	Gene vs Gene
		perl compare-genome.pl -WorkPath <path> -datatype <fragment> -seqtype <gene> -reference <reference gene> -isolate <isolate gene> -StorePath <output path>
	
	Protein vs Protein
		perl compare-genome.pl -WorkPath <path> -datatype <fragment> -seqtype <protein> -reference <reference protein> -isolate <isolate protein> -StorePath <output path>
	
	Genome vs Sequencing data (without amino acid annotation)
		perl compare-genome.pl -WorkPath <path> -datatype <sequencing> -reference <reference genome> -FQ1 <FASTQ R1> -FQ2 <FASTQ R2> -StorePath <output path>
		
	Genome vs Sequencing data (with amino acid annotation)
		perl compare-genome.pl -WorkPath <path> -datatype <sequencing> -reference <reference genome> -FQ1 <FASTQ R1> -FQ2 <FASTQ R2> -gff <gff file> -codon <codon table> -StorePath <output path>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-datatype	[str]	Input the data type (fragment or sequencing)
	
	-seqtype	[str]	Input the sequence type (genome, gene or protein), if the datatype fragment was chosen
	
	-reference	[str]	Input the reference sequence
	
	-isolate	[str]	Input the isolate sequence, if the datatype fragment was chosen
	
	-FQ1	[str]	Input the FASTQ R1 file, if the datatype sequencing was chosen
	
	-FQ2	[str]	Input the FASTQ R2 file, if the datatype sequencing was chosen
	
	-gff	[str]	Input the gene annotation file
	
	-codon	[int]	Input the translation table
	
	-StorePath	[str] Input the output path
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.04.22 00:53 0.0.1
	
	2021.05.51 15:22 0.0.2
	
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $trimmomatic = "/home/dell/software/Trimmomatic/trimmomatic-0.36.jar";
my $adapter = "/home/dell/software/Trimmomatic/adapters/TruSeq3-PE.fa";
my $bowtie2_build = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2-build";
my $bowtie2 = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2";
my $samtools = "/home/dell/miniconda3/envs/wgs/bin/samtools";
my $bcftools = "/home/dell/miniconda3/envs/wgs/bin/bcftools";
my $mafft = "/home/dell/miniconda3/envs/wgs/bin/mafft";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $prepar_seq = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/prepar_seq.pl";
my $compariative = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/comparativegenome.pl";
my $ATCG_SNP = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/ATCG_SNP.pl";
my $SNP2AA = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/SNP2AAMut.pl";
my $mergemut = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/mergemutation.pl";
my $genecomp = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/genecomparasion.pl";
my $protcomp = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/protcomparasion.pl";
my $ngscomp = "/var/www/other_scripts/bio-pipelines/Compare-genome/Functions/ngscomparative.pl";


my $analysis_path = "";
my ($datatype, $seqtype, $reference, $isolate, $gff, $help) = ("", "", "", "");
my $codon_tbl = 1;
my ($Fq_R1, $Fq_R2) = ("", "");
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'datatype=s' => \$datatype,
	'seqtype=s' => \$seqtype,
	'reference=s' => \$reference,
	'isolate:s' => \$isolate,
	'FQ1:s' => \$Fq_R1,
	'FQ2:s' => \$Fq_R2,
	'gff:s' => \$gff,
	'codon:i' => \$codon_tbl,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$datatype) or (!$storage_path) or ($help));
die `pod2text $0` if(($datatype eq "fragment") and ((!$reference) or (!$isolate) or (!$seqtype)));
die `pod2text $0` if(($datatype eq "sequencing") and ((!$reference) or (!$Fq_R1) or (!$Fq_R2)));

my @shell;
push @shell, "#!/usr/bin/bash\n";
if($datatype eq "fragment"){
	push @shell, $perl." ".$prepar_seq." -WorkPath ".$analysis_path." -reference ".$reference." -isolate ".$isolate."\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $mafft." --thread 20 --quiet ".$analysis_path."compare-sequences.fasta > ".$analysis_path."compare-sequences_maf.fasta\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $perl." ".$compariative." -WorkPath ".$analysis_path." -refname Ref -aligned ".$analysis_path."compare-sequences_maf.fasta\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	if($seqtype ne "protein"){
		push @shell, $perl." ".$ATCG_SNP." ".$analysis_path."snp.tsv > ".$analysis_path."ATCG_snp.tsv\n";
	}
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	if($gff and $seqtype eq "genome"){
		if(!$codon_tbl){
			push @shell, $perl." ".$SNP2AA." -WorkPath ".$analysis_path." -reference ".$reference." -gff ".$gff." -snp ATCG_snp.tsv\n";
		}
		else{
			push @shell, $perl." ".$SNP2AA." -WorkPath ".$analysis_path." -reference ".$reference." -gff ".$gff." -snp ATCG_snp.tsv -codon ".$codon_tbl."\n";
		}
		push @shell, $perl." ".$mergemut." -WorkPath ".$analysis_path." -snp ATCG_snp.tsv -del deletion.tsv -ins insertion.tsv -AAchange SNPAA.tsv -gff ".$gff."\n";
	}
	elsif(!$gff and $seqtype eq "genome"){
		push @shell, $perl." ".$mergemut." -WorkPath ".$analysis_path." -snp ATCG_snp.tsv -del deletion.tsv -ins insertion.tsv\n";
	}
	elsif(!$gff and $seqtype eq "gene"){
		if($codon_tbl){
			push @shell, $perl." ".$genecomp." -WorkPath ".$analysis_path." -reference ".$reference." -snp ATCG_snp.tsv -del deletion.tsv -ins insertion.tsv\n";
		}
		else{
			push @shell, $perl." ".$genecomp." -WorkPath ".$analysis_path." -reference ".$reference." -snp ATCG_snp.tsv -del deletion.tsv -ins insertion.tsv -codon ".$codon_tbl."\n";
		}
	}
	elsif(!$gff and $seqtype eq "protein"){
		push @shell, $perl." ".$protcomp." -WorkPath ".$analysis_path." -snp snp.tsv -del deletion.tsv -ins insertion.tsv\n";
	}
	push @shell, "chmod -R 777 ".$analysis_path."\n";
}
elsif($datatype eq "sequencing"){
	my $name = (split/_R1\./,$Fq_R1)[0];
	push @shell, "mkdir ".$analysis_path."Ref && chmod 777 -R ".$analysis_path."Ref/ && cp ".$analysis_path.$reference." ".$analysis_path."Ref\n";
	push @shell, "chmod 777 -R ".$analysis_path."Ref/"."\n";
	push @shell, $bowtie2_build." ".$analysis_path."Ref/".$reference." --thread 15 ".$analysis_path."Ref/Ref_bt2idx\n";
	push @shell, "chmod 777 -R ".$analysis_path."Ref/\n";
	push @shell, "mkdir ".$analysis_path."FASTQ && mkdir ".$analysis_path."Cleandata && mkdir ".$analysis_path."Trimdata\n";
	push @shell, "chmod 777 -R ".$analysis_path."FASTQ && chmod 777 -R ".$analysis_path."Cleandata && chmod 777 -R ".$analysis_path."Trimdata\n";
	push @shell, "mv ".$analysis_path.$Fq_R1." ".$analysis_path."FASTQ/"."\n";
	push @shell, "mv ".$analysis_path.$Fq_R2." ".$analysis_path."FASTQ/"."\n";
	push @shell, "chmod -R 777 ".$analysis_path."FASTQ"."\n";
	push @shell, "java -jar ".$trimmomatic." PE -threads 12 ".$analysis_path."FASTQ/".$Fq_R1." ".$analysis_path."FASTQ/".$Fq_R2." ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R1.fastq.gz ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R2.fastq.gz ILLUMINACLIP:".$adapter.":2:20:10:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:31\n";
	push @shell, "chmod -R 777 ".$analysis_path."Cleandata && chmod 777 -R ".$analysis_path."Trimdata\n";
	push @shell, "mkdir ".$analysis_path."Mapping\n";
	push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
	push @shell, $bowtie2." -p 30 -x ".$analysis_path."Ref/Ref_bt2idx -1 ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz -2 ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz -S ".$analysis_path."Mapping/".$name.".sam &> ".$analysis_path."Mapping/bt2.log\n";
	push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
	push @shell, $samtools." sort -@ 30 -o ".$analysis_path."Mapping/".$name.".bam ".$analysis_path."Mapping/".$name.".sam\n";
	push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
	push @shell, "rm -rf ".$analysis_path."Mapping/".$name.".sam\n";
	push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
	push @shell, $samtools." view -bF 4 -@ 30 ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name."_aligned.bam\n";
	push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
	push @shell, $samtools." index ".$analysis_path."Mapping/".$name."_aligned.bam > ".$analysis_path."Mapping/".$name."_aligned.bam.bai\n";
	push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
	push @shell, "mkdir ".$analysis_path."SNPcalls\n";
	push @shell, "chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."SNPcalls\n";
	push @shell, $samtools." mpileup -d 250 -gsf ".$analysis_path."Ref/".$reference." ".$analysis_path."Mapping/".$name."_aligned.bam -o ".$analysis_path."SNPcalls/".$name.".bcf\n";
	push @shell, "chmod 777 -R ".$analysis_path."SNPcalls\n";
	push @shell, $bcftools." view ".$analysis_path."SNPcalls/".$name.".bcf > ".$analysis_path."SNPcalls/".$name.".vcf\n";
	push @shell, "chmod 777 -R ".$analysis_path."SNPcalls\n";
	push @shell, $bcftools." call -mv ".$analysis_path."SNPcalls/".$name.".vcf > ".$analysis_path."SNPcalls/".$name."_SNPcalls.vcf\n";
	push @shell, "mkdir ".$analysis_path."Compargenome\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "cp ".$analysis_path.$reference." ".$analysis_path."Compargenome\n";
	push @shell, "cp ".$analysis_path."SNPcalls/".$name."_SNPcalls.vcf ".$analysis_path."Compargenome\n";
	push @shell, "chmod -R 777 ".$analysis_path."Compargenome\n";
	if($gff){
		push @shell, "cp ".$analysis_path.$gff." ".$analysis_path."Compargenome\n";
		push @shell, "chmod -R 777 ".$analysis_path."Compargenome\n";
		push @shell, $perl." ".$ngscomp." -WorkPath ".$analysis_path."Compargenome/ -reference ".$reference." -vcf ".$name."_SNPcalls.vcf -gff ".$gff." -codon ".$codon_tbl."\n";
		push @shell, "chmod -R 777 ".$analysis_path."Compargenome\n";
	}
	else{
		push @shell, $perl." ".$ngscomp." -WorkPath ".$analysis_path."Compargenome/ -reference ".$reference." -vcf ".$name."_SNPcalls.vcf\n";
		push @shell, "chmod -R 777 ".$analysis_path."Compargenome\n";
	}
	push @shell, "cp ".$analysis_path."Compargenome/comparative-results.tsv ".$analysis_path."\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
}
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."comparative-results.tsv ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
