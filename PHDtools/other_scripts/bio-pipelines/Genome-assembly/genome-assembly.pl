#!/usr/bin/env/perl -w

=head1 Description
	Target pathogen genome assemble
Usage
	perl genome-assembly.pl -WorkPath <analysis path> -FQ1 <FASTQ R1> -FQ2 <FASTQ R2> -Ref <reference genome> -StorePath <output path>
Parameters
	-WorkPath	[str]	Input the analysis path
	-FQ1	[str]	Input the FASTQ R1 file
	-FQ2	[str]	Input the FASTQ R2 file
	-Ref	[str]	Input the reference genome
	-taxid	[int]	Input the species level taxonomy ID of the target pathogen, optional, default = no
	-StorePath	[str] Input the output path
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2021.08.25 13:56 0.0.1
=cut

use strict;
use warnings;
use Getopt::Long;

my $ntdb = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/nt/nt";
my $acc2taxid_db = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/taxidb/acc_taxid_map";

my $trimmomatic = "/home/dell/software/Trimmomatic/trimmomatic-0.36.jar";
my $adapter = "/home/dell/software/Trimmomatic/adapters/TruSeq3-PE.fa";
my $bam2fastq = "/home/dell/software/bam2fastq-1.1.0/bam2fastq";
my $bowtie2_build = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2-build";
my $bowtie2 = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2";
my $samtools = "/home/dell/miniconda3/envs/wgs/bin/samtools";
my $megahit = "/home/dell/miniconda3/envs/wgs/bin/megahit";
my $makeblastdb = "/home/dell/software/ncbi-blast/bin/makeblastdb";
my $blastn = "/home/dell/software/ncbi-blast/bin/blastn";
my $perl= "/home/dell/miniconda3/envs/Perl/bin/perl";
my $extract_seq_pl = "/var/www/other_scripts/bio-pipelines/Genome-assembly/Functions/extractseq.pl";
my $sort_seq_pl = "/var/www/other_scripts/bio-pipelines/Genome-assembly/Functions/sortcontigs.pl";
my $tidy_ctg = "/var/www/other_scripts/bio-pipelines/Genome-assembly/Functions/tidyctgs.pl";
my $acc2taxid = "/var/www/other_scripts/bio-pipelines/Genome-assembly/Functions/Acc2taxid.pl";
my $extrac_taxidseq = "/var/www/other_scripts/bio-pipelines/Genome-assembly/Functions/extrataxseq.pl";

my $analysis_path = "";
my ($Fq_R1, $Fq_R2, $Ref, $help) = ("", "", "", "");
my $storage_path = "";
my $taxid = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'FQ1=s' => \$Fq_R1,
	'FQ2=s' => \$Fq_R2,
	'Ref=s' => \$Ref,
	'taxid:i' => \$taxid,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$Fq_R1) or (!$Fq_R2) or (!$Ref) or (!$storage_path) or ($help));
my $name = (split/_R1\./,$Fq_R1)[0];

my @shell;
push @shell, "#!/usr/bin/bash\n";
push @shell, "mkdir ".$analysis_path."Ref && chmod 777 ".$analysis_path."Ref/ && cp ".$analysis_path.$Ref." ".$analysis_path."Ref\n";
push @shell, "chmod 777 ".$analysis_path."Ref/".$Ref."\n";
push @shell, $bowtie2_build." ".$analysis_path."Ref/".$Ref." --thread 15 ".$analysis_path."Ref/Ref_bt2idx\n";
push @shell, $makeblastdb." -in ".$analysis_path."Ref/".$Ref." -dbtype nucl\n";
push @shell, "chmod 777 ".$analysis_path."Ref/*\n";
push @shell, "mkdir ".$analysis_path."FASTQ && mkdir ".$analysis_path."Cleandata && mkdir ".$analysis_path."Trimdata\n";
push @shell, "chmod 777 ".$analysis_path."FASTQ && chmod 777 ".$analysis_path."Cleandata && chmod 777 ".$analysis_path."Trimdata\n";
push @shell, "mv ".$analysis_path."*gz ".$analysis_path."FASTQ\n";
push @shell, "chmod -R 777 ".$analysis_path."FASTQ/\n";
push @shell, "java -jar ".$trimmomatic." PE -threads 12 ".$analysis_path."FASTQ/".$Fq_R1." ".$analysis_path."FASTQ/".$Fq_R2." ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R1.fastq.gz ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R2.fastq.gz ILLUMINACLIP:".$adapter.":2:20:10:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51\n";
push @shell, "chmod -R 777 ".$analysis_path."Cleandata && chmod 777 -R ".$analysis_path."Trimdata\n";
push @shell, "mkdir ".$analysis_path."Mapping && mkdir ".$analysis_path."Aligned_FQ && mkdir ".$analysis_path."Assemble\n";
push @shell, "chmod 777 ".$analysis_path."Mapping && chmod 777 ".$analysis_path."Aligned_FQ && chmod 777 ".$analysis_path."Assemble\n";
push @shell, $bowtie2." -p 30 -x ".$analysis_path."Ref/Ref_bt2idx -1 ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz -2 ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz -S ".$analysis_path."Mapping/".$name.".sam &> ".$analysis_path."Mapping/bt2.log\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." sort -@ 30 -o ".$analysis_path."Mapping/".$name.".bam ".$analysis_path."Mapping/".$name.".sam\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, "rm -rf ".$analysis_path."Mapping/".$name.".sam\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." flagstat ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name.".flagstat\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $samtools." view -bF 4 -@ 30 ".$analysis_path."Mapping/".$name.".bam > ".$analysis_path."Mapping/".$name."_aligned.bam\n";
push @shell, "chmod 777 -R ".$analysis_path."Mapping/\n";
push @shell, $bam2fastq." --aligned ".$analysis_path."Mapping/".$name."_aligned.bam -o ".$analysis_path."Aligned_FQ/".$name."_aligned#.fastq\n";
push @shell, "chmod 777 -R ".$analysis_path."Aligned_FQ/\n";
push @shell, "cat ".$analysis_path."Aligned_FQ/".$name."_aligned_1.fastq | gzip - > ".$analysis_path."Aligned_FQ/".$name."_aligned_R1.fastq.gz\n";
push @shell, "cat ".$analysis_path."Aligned_FQ/".$name."_aligned_2.fastq | gzip - > ".$analysis_path."Aligned_FQ/".$name."_aligned_R2.fastq.gz\n";
push @shell, "chmod -R 777 ".$analysis_path."Aligned_FQ/\n";
push @shell, $megahit." -m 32 -t 20 -1 ".$analysis_path."Aligned_FQ/".$name."_aligned_R1.fastq.gz -2 ".$analysis_path."Aligned_FQ/".$name."_aligned_R2.fastq.gz -o ".$analysis_path."Assemble/".$name." &> ".$analysis_path."Assemble/megahit.log\n";
push @shell, "chmod -R 777 ".$analysis_path."Assemble/ && chmod -R 777 ".$analysis_path."Assemble/".$name."\n";
push @shell, $blastn." -query ".$analysis_path."Assemble/".$name."/final.contigs.fa -db ".$analysis_path."Ref/".$Ref." -max_hsps 1 -num_threads 15 -evalue 1e-7 -outfmt 6 -out ".$analysis_path."Assemble/".$name."_blast.m6\n";
push @shell, "chmod -R 777 ".$analysis_path."Assemble/\n";
push @shell, "awk '{if(\$4>100)print\$1}' ".$analysis_path."Assemble/".$name."_blast.m6 | sort | uniq > ".$analysis_path."Assemble/".$name."_megahit_ID.txt\n";
push @shell, "chmod -R 777 ".$analysis_path."Assemble/\n";
push @shell, $perl." ".$extract_seq_pl." ".$analysis_path."Assemble/".$name."/final.contigs.fa ".$analysis_path."Assemble/".$name."_megahit_ID.txt > ".$analysis_path."Assemble/".$name."_extracted.fasta\n";
push @shell, "chmod -R 777 ".$analysis_path."Assemble/\n";
push @shell, $perl." ".$sort_seq_pl." megahit ".$analysis_path."Assemble/".$name."_extracted.fasta > ".$analysis_path."Assemble/".$name."_megahit_draftgenome.fasta\n";
push @shell, "chmod 777 -R ".$analysis_path."Assemble/\n";
push @shell, "cp ".$analysis_path."Assemble/".$name."_megahit_draftgenome.fasta ".$analysis_path."Final_Assembly.fasta\n";
push @shell, $perl." ".$tidy_ctg." -WorkPath ".$analysis_path." -contigs Final_Assembly.fasta"."\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "rm -rf ".$analysis_path."Final_Assembly.fasta\n";
push @shell, "mv ".$analysis_path."Final_assembly.fasta ".$analysis_path."Final_Assembly.fasta"."\n";
if($taxid){
	push @shell, "mkdir ".$analysis_path."Annotation"."\n";
	push @shell, $blastn." -query ".$analysis_path."Final_Assembly.fasta -db ".$ntdb." -max_hsps 1 -max_target_seqs 1 -num_threads 15 -evalue 1e-7 -outfmt 6 -out ".$analysis_path."Annotation/".$name."_blastn\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $perl." ".$acc2taxid." ".$analysis_path."Annotation/".$name."_blastn ".$acc2taxid_db."\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "mv ".$analysis_path."Annotation/".$name."_acc2taxid.map ".$analysis_path."\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $perl." ".$extrac_taxidseq." -WorkPath ".$analysis_path." -taxidtbl ".$name."_acc2taxid.map -taxid ".$taxid." -assemble Final_Assembly.fasta\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "rm -rf ".$analysis_path."Final_Assembly.fasta\n";
	push @shell, "mv ".$analysis_path."Final_assembly.fasta ".$analysis_path."Final_Assembly.fasta"."\n";
}
push @shell, "chmod 777 -R ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."Final_Assembly.fasta ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
