#!/usr/bin/perl -w

=head1 Description
	Report the quality of the whole genome sequencing result of the target pathogen
Usage
	perl wgs-quality.pl -WorkPath <path> -FQ1 <R1> -FQ2 <R2> -REF <refseq> -StorePath <output path>
Parameters
	-WorkPath		[str]	Analysis path, required
	-FQ1		[str]	FASTQ R1, required
	-FQ2		[str]	FASTQ R2, required
	-REF		[str]	Reference genome, required
	-StorePath	[str] Input the output path
	-h/-help		[str]	print help
Author
	Dongyan Xiong
Edit Time
	2022.07.14 16:01 0.0.1
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
my $bedtools = "/home/dell/miniconda3/envs/wgs/bin/bedtools";
my $bam2fastq = "/home/dell/software/bam2fastq-1.1.0/bam2fastq";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $fqstat = "/var/www/other_scripts/bio-pipelines/Wgs-quality/Functions/fqstat.pl";
my $covANDdep = "/var/www/other_scripts/bio-pipelines/Wgs-quality/Functions/covANDdep.pl";

my $analysis_path = "";
my ($Fq_R1, $Fq_R2, $Ref) = ("", "", "");
my $storage_path = "";
my $help = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'FQ1=s' => \$Fq_R1,
	'FQ2=s' => \$Fq_R2,
	'REF=s' => \$Ref,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);

my $name = (split/_R1\./,$Fq_R1)[0];
die `pod2text $0` if((!$analysis_path) or (!$Fq_R1) or (!$Fq_R2) or (!$Ref) or (!$storage_path) or ($help));

my @shell;
push @shell, "#!/usr/bin/bash\n";
push @shell, "mkdir ".$analysis_path."Ref && chmod 777 ".$analysis_path."Ref/ && cp ".$analysis_path.$Ref." ".$analysis_path."Ref\n";
push @shell, "chmod 777 ".$analysis_path."Ref/".$Ref."\n";
push @shell, $bowtie2_build." ".$analysis_path."Ref/".$Ref." --thread 15 ".$analysis_path."Ref/Ref_bt2idx\n";
push @shell, "chmod 777 ".$analysis_path."Ref/*\n";
push @shell, "mkdir ".$analysis_path."FASTQ && mkdir ".$analysis_path."Cleandata && mkdir ".$analysis_path."Trimdata\n";
push @shell, "mkdir ".$analysis_path."Quality/\n";
push @shell, "chmod -R 777 ".$analysis_path."Quality/\n";
push @shell, "chmod 777 ".$analysis_path."FASTQ && chmod 777 ".$analysis_path."Cleandata && chmod 777 ".$analysis_path."Trimdata\n";
push @shell, "mv ".$analysis_path."*gz ".$analysis_path."FASTQ\n";
push @shell, "chmod -R 777 ".$analysis_path."FASTQ/\n";
push @shell, $perl." ".$fqstat." -FQ1 ".$analysis_path."FASTQ/".$Fq_R1." -FQ2 ".$analysis_path."FASTQ/".$Fq_R2." > ".$analysis_path."Quality/".$name."-raw-reads.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Quality/\n";
push @shell, "java -jar ".$trimmomatic." PE -threads 12 ".$analysis_path."FASTQ/".$Fq_R1." ".$analysis_path."FASTQ/".$Fq_R2." ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R1.fastq.gz ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R2.fastq.gz ILLUMINACLIP:".$adapter.":2:20:10:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51\n";
push @shell, "chmod -R 777 ".$analysis_path."Cleandata && chmod 777 -R ".$analysis_path."Trimdata\n";
push @shell, $perl." ".$fqstat." -FQ1 ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz -FQ2 ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz > ".$analysis_path."Quality/".$name."-clean-reads.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Quality/\n";
push @shell, "mkdir ".$analysis_path."Mapping && mkdir ".$analysis_path."Aligned_FQ\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, $bowtie2." -p 30 -x ".$analysis_path."Ref/Ref_bt2idx -1 ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz -2 ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz -S ".$analysis_path."Mapping/".$name.".sam &> ".$analysis_path."Mapping/bt2.log\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." sort -@ 30 -o ".$analysis_path."Mapping/".$name.".bam ".$analysis_path."Mapping/".$name.".sam\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." flagstat ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name.".flagstat\n";
push @shell, $samtools." view -bF 4 -@ 30 ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name."_aligned.bam\n";
push @shell, "rm -rf ".$analysis_path."Mapping/".$name.".sam\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." index ".$analysis_path."Mapping/".$name."_aligned.bam > ".$analysis_path."Mapping/".$name."_aligned.bam.bai\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $bedtools." genomecov -ibam ".$analysis_path."Mapping/".$name."_aligned.bam -d > ".$analysis_path."Mapping/".$name.".beddep\n";
push @shell, $perl." ".$covANDdep." -genome ".$analysis_path."Ref/".$Ref." -DEPfile ".$analysis_path."Mapping/".$name.".beddep > ".$analysis_path."Quality/".$name."_covANDdep.tsv\n";
push @shell, "mkdir ".$analysis_path."Aligned_FQ"."\n";
push @shell, "chmod 777 -R ".$analysis_path."Aligned_FQ/\n";
push @shell, $bam2fastq." --aligned ".$analysis_path."Mapping/".$name."_aligned.bam -o ".$analysis_path."Aligned_FQ/".$name."_aligned#.fastq\n";
push @shell, "chmod 777 -R ".$analysis_path."Aligned_FQ/\n";
push @shell, $perl." ".$fqstat." -FQ1 ".$analysis_path."Aligned_FQ/".$name."_aligned_1.fastq -FQ2 ".$analysis_path."Aligned_FQ/".$name."_aligned_2.fastq > ".$analysis_path."Quality/".$name."-target-pair-reads.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Quality/\n";
push @shell, "zip -qjr ".$analysis_path."Quality_Reports.zip ".$analysis_path."Quality/".$name."-raw-reads.tsv ".$analysis_path."Quality/".$name."-clean-reads.tsv ".$analysis_path."Quality/".$name."-target-pair-reads.tsv ".$analysis_path."Quality/".$name."_covANDdep.tsv\n";
push @shell, "chmod 777 -R ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."Quality_Reports.zip ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
