#!/usr/bin/perl -w

=head1 Description
	Strain typing
Usage
	perl strain-typing.pl -WorkPath <analysis path> -model <AlleleCalling or WholeGenome> -StorePath <output path>
Parameters
	-WorkPath	[str]	Input the analysis path
	-model	[str]	AlleleCalling or WholeGenome
	-len	[int]	Input the gene length in AlleleCalling model, default = 200
	-percent	[float]	Input the gene present ratio in a population in AlleleCalling model, default = 0.99
	-codon	[int]	Input the translation table in AlleleCalling model, default = 1
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

my $prodigal = "/home/dell/software/Prodigal/bin/prodigal";
my $makeblastdb = "/home/dell/software/ncbi-blast/bin/makeblastdb";
my $blastn = "/home/dell/software/ncbi-blast/bin/blastn";
my $mafft = "/home/dell/miniconda3/envs/wgs/bin/mafft";
my $fasttree = "/home/dell/software/FastTree/FastTree";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $mergeallele = "/var/www/other_scripts/bio-pipelines/Strain-typing/Functions/merge-panalleles.pl";
my $tidymaf = "/var/www/other_scripts/bio-pipelines/Strain-typing/Functions/tidy-mafseq.pl";
my $linkallele = "/var/www/other_scripts/bio-pipelines/Strain-typing/Functions/link-panallele.pl";
my $seqparpre = "/var/www/other_scripts/bio-pipelines/Strain-typing/Functions/seqprepar.pl";

my $analysis_path = "";
my ($model, $genelen, $percent, $transtable, $help) = ("", 200, 0.99, 1, "");
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'model=s' => \$model,
	'len:i' => \$genelen,
	'percent:f' => \$percent,
	'codon:i' => \$transtable,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$model) or (!$storage_path) or ($help));
die `pod2text $0` if($model ne "AlleleCalling" and $model ne "WholeGenome");

my @shell;
push @shell, "#!/usr/bin/bash\n";
if($model eq "AlleleCalling"){
	my ($refgenome, $genomes) = ($analysis_path."refgenome", $analysis_path."genomes.fa");
	push @shell, "mkdir ".$analysis_path."genes_predict && mkdir ".$analysis_path."genes_gff && mkdir ".$analysis_path."all-to-all-blast && mkdir ".$analysis_path."merge-alleles\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "ls ".$analysis_path."genomes.fa | while read -r line; do ".$prodigal." -i ".$analysis_path."genomes.fa/\${line} -g ".$transtable." -d ".$analysis_path."genes_predict/\${line%.*}.fna -f gff -o ".$analysis_path."genes_gff/\${line%.*}.gff; done\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."genes_predict && chmod -R 777 ".$analysis_path."genes_gff\n";
	push @shell, "ls ".$analysis_path."genes_predict | while read -r line; do ".$makeblastdb." -in ".$analysis_path."genes_predict/\${line} -dbtype nucl; done\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."genes_predict\n";
	push @shell, "du -ha ".$analysis_path."genes_predict/*fna | sort -rn | head -n 1 | cut -f 2 | xargs -i cp {} ".$analysis_path."refgenome/reference.fna\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."refgenome\n";
	push @shell, "ls ".$analysis_path."genomes.fa | while read -r line; do ".$blastn." -query ".$analysis_path."refgenome/reference.fna -db ".$analysis_path."genes_predict/\${line%.*}.fna -max_target_seqs 10 -max_hsps 1 -evalue 1e-8 -num_threads 15 -outfmt 6 -out ".$analysis_path."all-to-all-blast/\${line%.*}.m6; done\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."all-to-all-blast\n";
	push @shell, $perl." ".$mergeallele." -reference ".$analysis_path."refgenome/reference.fna -minLen ".$genelen." -isolates ".$analysis_path."genes_predict -blastres ".$analysis_path."all-to-all-blast -presentrate ".$percent." -outpath ".$analysis_path."merge-alleles\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."merge-alleles\n";
	push @shell, "ls ".$analysis_path."merge-alleles/*sets.fasta | while read -r line; do ".$mafft." --thread 20 --quiet \${line} > \${line%.*}_maf.fasta; done\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."merge-alleles\n";
	push @shell, "ls ".$analysis_path."merge-alleles/*sets_maf.fasta | while read -r line; do ".$perl." ".$tidymaf." \${line} \${line%.*}-tidy.fasta; done\n";
	push @shell, "chmod -R 777 ".$analysis_path." && chmod -R 777 ".$analysis_path."merge-alleles\n";
	my $analysis_path2 = $analysis_path;
	chop($analysis_path2);
	push @shell, $perl." ".$linkallele." ".$analysis_path."merge-alleles ".$analysis_path."genes_predict ".$analysis_path2."\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $fasttree." -nt -gtr -fastest ".$analysis_path."Total_link.fasta > ".$analysis_path."StrainTyping.tree\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
}
elsif($model eq "WholeGenome"){
	my $genomes = $analysis_path."genomes.fa";
	push @shell, "cat ".$genomes."/*fa* >> ".$analysis_path."Total_genomes.fasta\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $perl." ".$seqparpre." ".$analysis_path."Total_genomes.fasta > ".$analysis_path."Prepared_genomes.fasta\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $mafft." --thread 12 --quiet ".$analysis_path."Prepared_genomes.fasta > ".$analysis_path."Prepared_genomes_maf.fasta\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $fasttree." -nt -gtr -fastest ".$analysis_path."Prepared_genomes_maf.fasta > ".$analysis_path."StrainTyping.tree\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, "sed -i ".'"'."s/\.fasta//g".'" '.$analysis_path."StrainTyping.tree\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
}
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."StrainTyping.tree ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
