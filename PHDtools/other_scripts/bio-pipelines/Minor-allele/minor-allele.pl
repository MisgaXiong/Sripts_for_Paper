#!/usr/bin/env/perl -w

=head1 Description
	intra-host variation identification
=head1 Usage
	perl minor-allele.pl -WorkPath <path> -FQ1 <R1> -FQ2 <R2> -REF <refseq> -MINDEP <min depth> -MAF <minor allele freq> -FLAG <flag> -StorePath <output path>
=head1 Parameters
	-WorkPath		[str]	Analysis path, required
	-FQ1		[str]	FASTQ R1, required
	-FQ2		[str]	FASTQ R2, required
	-REF		[str]	Reference genome, required
	-MINDEP		[int]	Minimum number of sequencing depth, default = 100, recommoned >= 100, optional
	-MAF		[float]	Minimum number of minor allele frequency, default = 0.02, optional
	-FLAG		[int]	flag value, required
	-Q		[int]	reads aligned quality, optional
	-GFF		[str]	gene annotation file, if the amino acid change want to be calculated, optional
	-CODON		[int]	Translation table, if the -GFF was chosen, default = 1, optional
	-StorePath	[str] Input the output path
	-h/-help		[str]	print help
=head1 Author
	Dongyan Xiong
=head1 Edit Time
	2022.06.30 16:07 0.0.1
=cut

use strict;
use warnings;
use Getopt::Long;
 
my $trimmomatic = "/home/dell/software/Trimmomatic/trimmomatic-0.36.jar";
my $adapter = "/home/dell/software/Trimmomatic/adapters/TruSeq3-PE.fa";
my $bowtie2_build = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2-build";
my $bowtie2 = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2";
my $samtools = "/home/dell/miniconda3/envs/wgs/bin/samtools";
my $bcftools = "/home/dell/miniconda3/envs/wgs/bin/bcftools";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $alignstat = "/var/www/other_scripts/bio-pipelines/Minor-allele/Functions/AlignStatplus2.pl";
my $minorallele = "/var/www/other_scripts/bio-pipelines/Minor-allele/Functions/MinorAlleleCalling.pl";
my $snpfind = "/var/www/other_scripts/bio-pipelines/Minor-allele/Functions/SNPFind.pl";

my $analysis_path = "";
my ($Fq_R1, $Fq_R2, $Ref, $depth, $iSNV_rate) = ("", "", "", 100, 0.02);
my ($flag, $quality, $gff, $codon_tbl) = ("", "", "", 1);
my $storage_path = "";
my $help;
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'FQ1=s' => \$Fq_R1,
	'FQ2=s' => \$Fq_R2,
	'REF=s' => \$Ref,
	'MINDEP=i' => \$depth,
	'MAF=f' => \$iSNV_rate,
	'FLAG=i' => \$flag,
	'Q:i' => \$quality,
	'GFF:s' => \$gff,
	'CODON:i' => \$codon_tbl,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
my $name = (split/_R1\./,$Fq_R1)[0];
die `pod2text $0` if((!$analysis_path) or (!$Fq_R1) or (!$Fq_R2) or (!$Ref) or (!$depth) or (!$iSNV_rate) or (!$flag) or (!$storage_path) or ($help));

my @shell;
push @shell, "#!/usr/bin/bash\n";
push @shell, "mkdir ".$analysis_path."Ref && chmod 777 -R ".$analysis_path."Ref/ && cp ".$analysis_path.$Ref." ".$analysis_path."Ref\n";
push @shell, "chmod 777 -R ".$analysis_path."Ref/"."\n";
push @shell, $bowtie2_build." ".$analysis_path."Ref/".$Ref." --thread 15 ".$analysis_path."Ref/Ref_bt2idx\n";
push @shell, "chmod 777 ".$analysis_path."Ref/*\n";
if($gff){
	push @shell, "cp ".$analysis_path.$gff." ".$analysis_path."Ref\n";
	push @shell, "chmod 777 -R ".$analysis_path."Ref/"."\n";
}
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
push @shell, $samtools." flagstat ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name.".flagstat\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
if(!$quality){
	push @shell, $samtools." view -bF ".$flag." -@ 30 ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name."_aligned.bam\n";
}
else{
	push @shell, $samtools." view -bF ".$flag." -q ".$quality." -@ 30 ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name."_aligned.bam\n";
}
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." index ".$analysis_path."Mapping/".$name."_aligned.bam > ".$analysis_path."Mapping/".$name."_aligned.bam.bai\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, "mkdir ".$analysis_path."SNPcalls\n";
push @shell, "chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."SNPcalls\n";
push @shell, $samtools." mpileup -d 250 -gsf ".$analysis_path."Ref/".$Ref." ".$analysis_path."Mapping/".$name."_aligned.bam -o ".$analysis_path."SNPcalls/".$name.".bcf\n";
push @shell, "chmod 777 -R ".$analysis_path."SNPcalls\n";
push @shell, $bcftools." view ".$analysis_path."SNPcalls/".$name.".bcf > ".$analysis_path."SNPcalls/".$name.".vcf\n";
push @shell, "chmod 777 -R ".$analysis_path."SNPcalls\n";
push @shell, $bcftools." call -mv ".$analysis_path."SNPcalls/".$name.".vcf > ".$analysis_path."SNPcalls/".$name."_SNPcalls.vcf\n";
push @shell, "mkdir ".$analysis_path."BaseCount\n";
push @shell, "chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."BaseCount\n";
push @shell, $samtools." view ".$analysis_path."Mapping/".$name."_aligned.bam > ".$analysis_path."BaseCount/".$name."_bam.txt\n";
push @shell, $perl." ".$alignstat." -genome ".$analysis_path."Ref/".$Ref." -input ".$analysis_path."BaseCount/".$name."_bam.txt -output ".$analysis_path."BaseCount/".$name."-6base.tsv\n";
push @shell, "chmod 777 -R ".$analysis_path."BaseCount/\n";
push @shell, "rm -rf ".$analysis_path."BaseCount/".$name."_bam.txt\n";
push @shell, "chmod 777 -R ".$analysis_path."BaseCount/\n";
if(!$gff){
	push @shell, $perl." ".$minorallele." -reference ".$analysis_path."Ref/".$Ref." -baseProfile ".$analysis_path."BaseCount/".$name."-6base.tsv -mindepth ".$depth." -minMuF ".$iSNV_rate." -output ".$analysis_path."BaseCount/".$name."-iSNV-profile.tsv"."\n";
}
else{
	push @shell, $perl." ".$minorallele." -reference ".$analysis_path."Ref/".$Ref." -baseProfile ".$analysis_path."BaseCount/".$name."-6base.tsv -mindepth ".$depth." -minMuF ".$iSNV_rate." -gff ".$analysis_path."Ref/".$gff." -codon ".$codon_tbl." -output ".$analysis_path."BaseCount/".$name."-iSNV-profile.tsv"."\n";
}
push @shell, $perl." ".$snpfind." ".$analysis_path."Ref/".$Ref." ".$analysis_path."SNPcalls/".$name."_SNPcalls.vcf ".$analysis_path."BaseCount/".$name."-6base.tsv\n";
push @shell, "chmod 777 -R ".$analysis_path."BaseCount/\n";
push @shell, "zip -qjr ".$analysis_path."iSNV_Results.zip ".$analysis_path."BaseCount/".$name."-iSNV-profile.tsv ".$analysis_path."BaseCount/".$name."_poteintial_SNP.txt ".$analysis_path."BaseCount/".$name."_SNP.txt\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."iSNV_Results.zip ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
