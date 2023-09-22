#!/usr/bin/perl -w

=head1 Description
	PCR primers-probe aligment in the strain level
Usage
	perl primer-align.pl -WorkPath <path> -genomes <genome of target pathogen fasta> -primers <primers-probe fasta>
	
	or
	
	perl primer-align.pl -WorkPath <path> -genomes <genome of target pathogen fasta> -primers <primers-probe fasta> -amplicon <amplicon sequence fasta>
Parameters
	-WorkPath	[str]	Input the analysis path
	-genomes	[str]	Input the isolate genomes fasta
	-primers	[str]	Input the primers-probe fasta
	-amplicon	[str]	Input the amplicon sequence
	-StorePath	[str] Input the output path
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.18 22:56 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $makeblastdb = "/home/dell/software/ncbi-blast/bin/makeblastdb";
my $blastn = "/home/dell/software/ncbi-blast/bin/blastn";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $amplicon_find = "/var/www/other_scripts/bio-pipelines/Primer-align/Functions/primers2amplicon.pl";
my $extract_seq = "/var/www/other_scripts/bio-pipelines/Primer-align/Functions/extractseq2primeralign.pl";
my $primer_align = "/var/www/other_scripts/bio-pipelines/Primer-align/Functions/PrimerAlign.pl";

my $analysis_path = "";
my ($isolates, $PCR_primers, $amplicon, $help) = ("", "", "", "");
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'genomes=s' => \$isolates,
	'primers=s' => \$PCR_primers,
	'amplicon:s' => \$amplicon,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$isolates) or (!$PCR_primers) or (!$storage_path) or ($help));

my @shell;
push @shell, "#!/usr/bin/bash\n";
push @shell, $makeblastdb." -in ".$analysis_path.$isolates." -dbtype nucl -logfile ".$analysis_path."makedb.log\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
if(!$amplicon){
	push @shell, $blastn." -task blastn-short -query ".$analysis_path.$PCR_primers." -db ".$analysis_path.$isolates." -max_target_seqs 500000 -max_hsps 1 -num_threads 15 -outfmt 6 -out ".$analysis_path."rapid-align.tb6\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $perl." ".$amplicon_find." -WorkPath ".$analysis_path." -genomes ".$isolates." -primers ".$PCR_primers." -primerm6 rapid-align.tb6\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $blastn." -task dc-megablast -query ".$analysis_path."AmpliconSeq.fasta -db ".$analysis_path.$isolates." -max_target_seqs 500000 -max_hsps 1 -num_threads 15 -outfmt 6 -out ".$analysis_path."target-align.tb6\n";
}
else{
	push @shell, $blastn." -task dc-megablast -query ".$analysis_path.$amplicon." -db ".$analysis_path.$isolates." -max_target_seqs 500000 -max_hsps 1 -num_threads 15 -outfmt 6 -out ".$analysis_path."target-align.tb6\n";
	push @shell, $makeblastdb." -in ".$analysis_path.$amplicon." -dbtype nucl -logfile ".$analysis_path."makedb.log\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $blastn." -task blastn-short -query ".$analysis_path.$PCR_primers." -db ".$analysis_path.$amplicon." -max_target_seqs 500000 -max_hsps 1 -num_threads 15 -outfmt 6 -out ".$analysis_path."rapid-align.tb6\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
	push @shell, $perl." ".$amplicon_find." -WorkPath ".$analysis_path." -genomes ".$amplicon." -primers ".$PCR_primers." -primerm6 rapid-align.tb6\n";
	push @shell, "chmod -R 777 ".$analysis_path."\n";
}
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, $perl." ".$extract_seq." -WorkPath ".$analysis_path." -genomes ".$isolates." -targetm6 target-align.tb6\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, $perl." ".$primer_align." -WorkPath ".$analysis_path." -pmlist pmlist.tsv -target Sequences2PrimerAlign.fasta\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "zip -qjr ".$analysis_path."StrainLevel-PrimerAlignment.zip ".$analysis_path."CompletedMatch.tsv ".$analysis_path."UncompletedMatch.tsv ".$analysis_path."PrimerMutation-Matrix.tsv ".$analysis_path."RealPrimerMutation-Matrix.tsv ".$analysis_path."PrimerStat.tsv ".$analysis_path."lowquality-genomes-uncovered-Amplicon.txt\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."StrainLevel-PrimerAlignment.zip ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
