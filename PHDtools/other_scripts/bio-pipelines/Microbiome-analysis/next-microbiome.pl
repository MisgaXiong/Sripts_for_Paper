#!/usr/bin/env/perl -w

=head1 Description
	Pathogen X identification through mNGS approach
Usage
	perl next-microbiome.pl -WorkPath <analysis path> -FQ1 <FASTQ R1> -FQ2 <FASTQ R2> -Host <Host name> -BGM <Background microbiota sequences>
Parameters
	-WorkPath	[str]	Input the analysis path
	-FQ1	[str]	Input the FASTQ R1 file, if the datatype sequencing was chosen
	-FQ2	[str]	Input the FASTQ R2 file, if the datatype sequencing was chosen
	-Host	[str]	Input the Host name ("Human", "Pig", "UNKNOWN"), optional
	-BGM	[str]	Input the background microbiota genomes, optional
	-contigLen	[int]	Length of contigs to identify species
	-StorePath	[str]	Input the result storage path
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2021.05.30 14:56 0.0.1
	2021.08.22 11:34 0.0.2
=cut

use strict;
use warnings;
use Getopt::Long;

my %host_path = (
	"HUMAN" => "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/Genome_idx/Human/Genome_idx/Hsa_genome",
	"PIG" => "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/Genome_idx/Pig/Genome_idx/Ssc_genome",
	"MOUSE" => "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/Genome_idx/Mouse/Genome_idx/Mms_genome",
);
my $ntdb = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/nt/nt";
my $acc2taxid = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/taxidb/acc_taxid_map";
my $dmpnames = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/taxidb/taxdump/names.dmp";
my $TaxInfo = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/taxidb/taxdump/TaxInfo.txt";
my $gmelens = "/mnt/768e536c-582e-486a-9d82-95be22b89714/database/blastdb/taxidb/genomes.csv";

my $fastq_screen = "/home/dell/miniconda3/envs/fqscreen/share/fastq-screen-0.15.3-0/fastq_screen";
my $trimmomatic = "/home/dell/software/Trimmomatic/trimmomatic-0.36.jar";
my $adapter = "/home/dell/software/Trimmomatic/adapters/TruSeq3-PE.fa";
my $bowtie2 = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2";
my $bowtie2_build = "/home/dell/software/bowtie2-2.3.5-sra-linux-x86_64/bowtie2-build";
my $megahit = "/home/dell/miniconda3/envs/wgs/bin/megahit";
my $samtools = "/home/dell/miniconda3/envs/wgs/bin/samtools";
my $blastn = "/home/dell/software/ncbi-blast/bin/blastn";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $sortcontigs = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/sortcontigs.pl";
my $longercontigs = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/longerCongtigs.pl";
my $RPMvalue = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/RPMvalue.pl";
my $accTotaxid = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/Acc2taxid.pl";
my $taxid2name = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/Taxid2name.pl";
my $AnotSP = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/AnotSP.pl";
my $MetaProfile = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/MetaProfile.pl";
my $abundancesort = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/abundancesort.pl";
my $coverage = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/coverage.pl";
my $isolate = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/Functions/species-isolates.pl";

my $analysis_path = "";
my ($Fq_R1, $Fq_R2, $host, $bgm, $help) = ("", "", "UNKNOWN", "", "");
my $contig_len = 200;
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'FQ1=s' => \$Fq_R1,
	'FQ2=s' => \$Fq_R2,
	'Host:s' => \$host,
	'BGM:s' => \$bgm,
	'contigLen:i' => \$contig_len,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$Fq_R1) or (!$Fq_R2) or (!$storage_path) or ($help));
my $name = (split/_R1\./,$Fq_R1)[0];

my (@fqscreen, @shell);
push @shell, "#!/usr/bin/bash\n";
push @shell, "mkdir ".$analysis_path."FASTQ && chmod -R 777 ".$analysis_path."\n";
push @shell, "mv ".$analysis_path.$Fq_R1." ".$analysis_path."FASTQ/"."\n";
push @shell, "mv ".$analysis_path.$Fq_R2." ".$analysis_path."FASTQ/"."\n";
push @shell, "chmod -R 777 ".$analysis_path."FASTQ"."\n";
if((uc($host) eq "UNKNOWN" or uc($host) eq "NO") and !$bgm){
	push @shell, "mkdir ".$analysis_path."Cleandata && mkdir ".$analysis_path."Trimdata\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "java -jar ".$trimmomatic." PE -threads 12 ".$analysis_path."FASTQ/".$Fq_R1." ".$analysis_path."FASTQ/".$Fq_R2." ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R1.fastq.gz ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R2.fastq.gz ILLUMINACLIP:".$adapter.":2:20:10:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:31\n";
}
else{
	if(!$bgm){
		push @fqscreen, "BOWTIE2 ".$bowtie2."\n";
		push @fqscreen, "THREADS\t\t8\n";
		push @fqscreen, "DATABASE\t".uc($host)."\t".$host_path{uc($host)}."\n";
	}
	else{
		if(uc($host) ne "UNKNOWN" and uc($host) ne "NO" and $bgm){
			push @fqscreen, "BOWTIE2 ".$bowtie2."\n";
			push @fqscreen, "THREADS\t\t10\n";
			push @fqscreen, "DATABASE\t".uc($host)."\t".$host_path{uc($host)}."\n";
			push @shell, "mkdir ".$analysis_path."BGM\n";
			push @shell, "chmod -R 777 ".$analysis_path."\n";
			push @shell, "mv ".$analysis_path.$bgm." ".$analysis_path."BGM\n";
			push @shell, "chmod -R 777 ".$analysis_path."\n";
			push @shell, $bowtie2_build." ".$analysis_path."BGM/".$bgm." --thread 15 ".$analysis_path."BGM/BGM_bt2idx\n";
			push @shell, "chmod -R 777 ".$analysis_path."\n";
			push @fqscreen, "DATABASE\t"."BackGroundMicrobiota"."\t".$analysis_path."BGM/BGM_bt2idx"."\n";
		}
		elsif((uc($host) eq "UNKNOWN" or uc($host) eq "NO") and $bgm){
			push @shell, "mkdir ".$analysis_path."BGM\n";
			push @shell, "chmod -R 777 ".$analysis_path."\n";
			push @shell, "mv ".$analysis_path.$bgm." ".$analysis_path."BGM\n";
			push @shell, "chmod -R 777 ".$analysis_path."\n";
			push @shell, $bowtie2_build." ".$analysis_path."BGM/".$bgm." --thread 15 ".$analysis_path."BGM/BGM_bt2idx\n";
			push @shell, "chmod -R 777 ".$analysis_path."\n";
			push @fqscreen, "BOWTIE2 ".$bowtie2."\n";
			push @fqscreen, "THREADS\t\t10\n";
			push @fqscreen, "DATABASE\t"."BackGroundMicrobiota"."\t".$analysis_path."BGM/BGM_bt2idx"."\n";
		}
	}
	push @shell, "mkdir ".$analysis_path."Fliterdata && chmod -R 777 ".$analysis_path." && mkdir ".$analysis_path."FliterRecord && chmod -R 777 ".$analysis_path."\n";
	push @shell, $fastq_screen." --conf ".$analysis_path."fastq_screen.conf --aligner bowtie2 --threads 10 --nohits ".$analysis_path."FASTQ/".$Fq_R1." ".$analysis_path."FASTQ/".$Fq_R2." --outdir ".$analysis_path."FASTQ &> ".$analysis_path."FliterRecord/".$name."_fqscreen.log\n";
	push @shell, "chmod -R 777 ".$analysis_path."FASTQ"."\n";
	push @shell, "chmod -R 777 ".$analysis_path."FliterRecord"."\n";
	push @shell, "mv ".$analysis_path."FASTQ/".$name."_R1.tagged_filter.fastq.gz ".$analysis_path."Fliterdata"."\n";
	push @shell, "mv ".$analysis_path."FASTQ/".$name."_R2.tagged_filter.fastq.gz ".$analysis_path."Fliterdata"."\n";
	push @shell, "chmod -R 777 ".$analysis_path."Fliterdata"."\n";
	push @shell, "mv ".$analysis_path."FASTQ/*screen* ".$analysis_path."FliterRecord\n";
	push @shell, "rm -rf ".$analysis_path."FASTQ/*tagged*\n";
	push @shell, "chmod -R 777 ".$analysis_path."FliterRecord"."\n";
	push @shell, "mkdir ".$analysis_path."Cleandata && mkdir ".$analysis_path."Trimdata\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "java -jar ".$trimmomatic." PE -threads 12 ".$analysis_path."Fliterdata/".$name."_R1.tagged_filter.fastq.gz ".$analysis_path."Fliterdata/".$name."_R2.tagged_filter.fastq.gz ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R1.fastq.gz ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz ".$analysis_path."Trimdata/".$name."_unpaired_R2.fastq.gz ILLUMINACLIP:".$adapter.":2:20:10:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:31\n";
	push @shell, "chmod -R 777 ".$analysis_path."Cleandata && chmod 777 -R ".$analysis_path."Trimdata\n";
}
push @shell, "mkdir ".$analysis_path."Assembly && chmod -R 777 ".$analysis_path."Assembly\n";
push @shell, $megahit." -m 32 -t 20 --min-contig-len ".$contig_len." -1 ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz -2 ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz -o ".$analysis_path."Assembly/".$name." &> ".$analysis_path."Assembly/".$name."_megahit.log\n";
push @shell, "chmod -R 777 ".$analysis_path."Assembly && chmod -R 777 ".$analysis_path."Assembly/".$name."\n";
push @shell, $perl." ".$sortcontigs." ".$analysis_path."Assembly/".$name."/final.contigs.fa > ".$analysis_path."Assembly/".$name."/final_assembly.sort.fasta\n";
push @shell, "chmod -R 777 ".$analysis_path."Assembly/".$name."\n";
push @shell, $perl." ".$longercontigs." 0 ".$analysis_path."Assembly/".$name."/final_assembly.sort.fasta > ".$analysis_path."Assembly/".$name."/Final_Assembly.fasta\n";
push @shell, "chmod -R 777 ".$analysis_path."Assembly/".$name."\n";
push @shell, "mkdir ".$analysis_path."Mapping && chmod -R 777 ".$analysis_path."Mapping\n";
push @shell, "mkdir ".$analysis_path."Mapping/".$name." && chmod -R 777 ".$analysis_path."Mapping/".$name."\n";
push @shell, "mkdir ".$analysis_path."Mapping/".$name."/Ref && mkdir ".$analysis_path."Mapping/".$name."/BAM\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."\n";
push @shell, $bowtie2_build." ".$analysis_path."Assembly/".$name."/Final_Assembly.fasta --thread 15 ".$analysis_path."Mapping/".$name."/Ref/".$name."_bt2idx\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/Ref\n";
push @shell, $bowtie2." -p 30 -x ".$analysis_path."Mapping/".$name."/Ref/".$name."_bt2idx -1 ".$analysis_path."Cleandata/".$name."_paired_R1.fastq.gz -2 ".$analysis_path."Cleandata/".$name."_paired_R2.fastq.gz -S ".$analysis_path."Mapping/".$name."/BAM/".$name.".sam &> ".$analysis_path."Mapping/".$name."/BAM/".$name."_bt2.log\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/BAM\n";
push @shell, $samtools." sort -@ 30 -o ".$analysis_path."Mapping/".$name."/BAM/".$name.".sort.bam ".$analysis_path."Mapping/".$name."/BAM/".$name.".sam\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/BAM && rm -rf ".$analysis_path."Mapping/".$name."/BAM/".$name.".sam\n";
push @shell, $samtools." index ".$analysis_path."Mapping/".$name."/BAM/".$name.".sort.bam\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/BAM\n";
push @shell, $samtools." idxstats ".$analysis_path."Mapping/".$name."/BAM/".$name.".sort.bam > ".$analysis_path."Mapping/".$name."/BAM/".$name.".stat\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/BAM\n";
push @shell, $perl." ".$RPMvalue." ".$analysis_path."Mapping/".$name."/BAM > ".$analysis_path."Mapping/".$name."/BAM/".$name."_RPM.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/BAM\n";
push @shell, "mkdir ".$analysis_path."Annotation && chmod -R 777 ".$analysis_path."Annotation\n";
push @shell, "mkdir ".$analysis_path."Annotation/".$name."\n";
push @shell, "chmod -R 777 ".$analysis_path."Annotation/".$name."\n";
push @shell, $blastn." -query ".$analysis_path."Assembly/".$name."/Final_Assembly.fasta -db ".$ntdb." -num_threads 15 -max_target_seqs 1 -max_hsps 1 -evalue 1e-8 -outfmt 6 -out ".$analysis_path."Annotation/".$name."/".$name."_blastn.m6\n";
push @shell, "chmod -R 777 ".$analysis_path."Annotation/".$name."\n";
push @shell, $perl." ".$accTotaxid." ".$analysis_path."Annotation/".$name."/".$name."_blastn.m6 ".$acc2taxid."\n";
push @shell, "chmod -R 777 ".$analysis_path."Annotation/".$name."\n";
push @shell, $perl." ".$taxid2name." ".$dmpnames." ".$analysis_path."Annotation/".$name."/".$name."_acc2taxid.map\n";
push @shell, "chmod -R 777 ".$analysis_path."Annotation/".$name."\n";
push @shell, $perl." ".$AnotSP." ".$TaxInfo." ".$analysis_path."Annotation/".$name."/".$name."_taxid2name.map\n";
push @shell, "chmod -R 777 ".$analysis_path."Annotation/".$name."\n";
push @shell, "mkdir ".$analysis_path."Metagenomic_Profile && chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$analysis_path."Metagenomic_Profile/".$name." && chmod -R 777 ".$analysis_path."Metagenomic_Profile\n";
push @shell, $perl." ".$MetaProfile." ".$analysis_path."Annotation/".$name."/".$name."_species.ant ".$analysis_path."Mapping/".$name."/BAM/".$name."_RPM.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Mapping/".$name."/BAM\n";
push @shell, "mv ".$analysis_path."Mapping/".$name."/BAM/".$name."_metagenomic.profile ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome.profile\n";
push @shell, "mv ".$analysis_path."Mapping/".$name."/BAM/".$name."_contigs_unidentify.tsv ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_contigs_unidentified.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Metagenomic_Profile/".$name."\n";
push @shell, $perl." ".$abundancesort." ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome.profile ".$TaxInfo."\n";
push @shell, "chmod -R 777 ".$analysis_path."Metagenomic_Profile/".$name."\n";
push @shell, $perl." ".$coverage." -reflens ".$gmelens." -microbes ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Metagenomic_Profile/".$name."\n";
push @shell, "rm -rf ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv\n";
push @shell, "mv ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv2 ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."Metagenomic_Profile/".$name."\n";
push @shell, "mkdir ".$analysis_path."Isolates && chmod -R 777 ".$analysis_path."Isolates\n";
push @shell, "mkdir ".$analysis_path."Isolates/".$name." && chmod -R 777 ".$analysis_path."Isolates\n";
push @shell, $perl." ".$isolate." ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv ".$analysis_path."Annotation/".$name."/".$name."_species.ant ".$analysis_path."Assembly/".$name."/Final_Assembly.fasta ".$analysis_path."Isolates/".$name."\n";
push @shell, "chmod -R 777 ".$analysis_path."Isolates/".$name."\n";
if((uc($host) eq "UNKNOWN" or uc($host) eq "NO") and !$bgm){
	push @shell, "zip -qjr ".$analysis_path."Microbiota_Results.zip ".$analysis_path."Isolates ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv\n";
}
else{
	push @shell, "zip -qjr ".$analysis_path."Microbiota_Results.zip ".$analysis_path."Isolates ".$analysis_path."Metagenomic_Profile/".$name."/".$name."_microbiome_sort.tsv ".$analysis_path."FliterRecord/".$name."_R1_screen.txt ".$analysis_path."FliterRecord/".$name."_R2_screen.txt\n";
}
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."Microbiota_Results.zip ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
if(@fqscreen){
	open(OUT, ">", $analysis_path."fastq_screen.conf");
	print OUT @fqscreen;
	close OUT;
}
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
